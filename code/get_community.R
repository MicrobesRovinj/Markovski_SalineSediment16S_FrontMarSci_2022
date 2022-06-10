#################################################################################################################
# get_community.R
#
# A script to format mothur taxonomy data.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.silva.wang.tax.summary
# Produces: results/numerical/community.Rdata
#
#################################################################################################################

# Loading data with sequence abundances and excluding eukaryotic sequences
community <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.silva.wang.tax.summary")
eukaryota <- filter(community, str_detect(taxon, "^Eukaryota$"))$rankID
community <- community %>%
  filter(!str_detect(rankID, paste0("^", eukaryota))) %>%
  filter(taxon != "Root")

# Removing chloroplast and mitochondrial sequences and subtracting their number from higher taxonomic levels to which they belong
chloroplast <- filter(community, str_detect(taxon, "^Chloroplast$"))$rankID
mitochondria <- filter(community, str_detect(taxon, "^Mitochondria$"))$rankID
community <- community %>%
  mutate(across(5 : ncol(.), ~ case_when(
  rankID == str_extract(chloroplast, "(\\d+\\.){3}\\d+") ~ . - .[taxon == "Chloroplast"],
  rankID == str_extract(chloroplast, "(\\d+\\.){2}\\d+") ~ . - .[taxon == "Chloroplast"],
  rankID == str_extract(chloroplast, "(\\d+\\.){1}\\d+") ~ . - .[taxon == "Chloroplast"],
  TRUE ~ .))) %>%
  filter(!str_detect(taxon, "^Chloroplast")) %>%
  mutate(across(5 : ncol(.), ~ case_when(
    rankID == str_extract(mitochondria, "(\\d+\\.){4}\\d+") ~ . - .[taxon == "Mitochondria"],
    rankID == str_extract(mitochondria, "(\\d+\\.){3}\\d+") ~ . - .[taxon == "Mitochondria"],
    rankID == str_extract(mitochondria, "(\\d+\\.){2}\\d+") ~ . - .[taxon == "Mitochondria"],
    rankID == str_extract(mitochondria, "(\\d+\\.){1}\\d+") ~ . - .[taxon == "Mitochondria"],
    TRUE ~ .))) %>%
  filter(!str_detect(taxon, "^Mitochondria")) %>%
  # Removing negative and positive controls from the data
  select(-ATCC_2, -ATCC_3, -ATCC_4, -ATCC_5, -NC_3, -NC_4) %>%
  group_by(taxlevel) %>%
  mutate_at(5 : ncol(.), ~ . / sum(.) * 100) %>%
  ungroup()
  
# Saving
save(community, file = "results/numerical/community.Rdata")
  