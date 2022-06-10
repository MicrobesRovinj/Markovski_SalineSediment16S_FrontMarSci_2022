#################################################################################################################
# get_estimators_indices_metadata.R
#
# A script to calculate richness estimators and diversity indices.
# Dependencies: results/numerical/rarefied.Rdata
#               data/raw/metadata.csv
# Produces: results/numerical/estimators_indices_metadata.Rdata
#
#################################################################################################################

# Loading the rarefied community data
load("results/numerical/rarefied.Rdata")

# Copying sample labels to rows names (input for library vegan)
rarefied <- rarefied %>%
  column_to_rownames("Group")

# Calculating observed number of OTUs and estimators Chao1 and ACE
estimators <- rarefied %>%
  estimateR() %>%
  t() %>%
  as_tibble(.name_repair = "check_unique", rownames = NA) %>%
  rownames_to_column("Group")

# Calculating diversity indices
shannon <- rarefied %>%
  diversity(index = "shannon") %>%
  enframe(name = "Group", value = "shannon")
invsimpson <- rarefied %>%
  diversity(index = "invsimpson") %>%
  enframe(name = "Group", value = "invsimpson")

# Transforming Shannon entropy to effective number of OTUs
# (http://www.loujost.com/Statistics%20and%20Physics/Diversity%20and%20Similarity/EffectiveNumberOfSpecies.htm)
eshannon <- mutate(shannon, shannon = exp(shannon)) %>%
  rename(eshannon = shannon)

# Joining estimators and indices
estimators_indices <- inner_join(estimators, eshannon, by = "Group") %>%
  inner_join(., invsimpson, by = "Group")

# Loading metadata
metadata <- read_tsv("data/raw/metadata.csv")

# Joining metadata with estimators and indices
Sys.setlocale(locale = "en_GB.utf8")
estimators_indices_metadata <- inner_join(metadata, estimators_indices, by = c("ID" = "Group")) %>%
  mutate(date = as.Date(date, "%d.%m.%Y"))

# Saving calculated estimators and indices
save(estimators_indices_metadata, file = "results/numerical/estimators_indices_metadata.Rdata")
