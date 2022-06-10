#################################################################################################################
# get_rarefied.R
#
# A script to create a randomly rarefied community data.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared
# Produces: results/numerical/rarefied.Rdata
#
#################################################################################################################

# Loading OTU/sample data
shared <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared")

# Adding sample labels to row names (input for library vegan)
shared <- shared %>%
  select(-label, -numOtus) %>%
  column_to_rownames("Group")

# Generating randomly rarefied community data
rarefied <- shared %>%
  rrarefy(., min(rowSums(.))) %>%
  as_tibble(.name_repair = "check_unique", rownames = NA) %>%
  rownames_to_column("Group") %>%
  select(!where(is.numeric) | where(~ is.numeric(.) && sum(.) != 0))

# Saving rarefied community data
save(rarefied, file = "results/numerical/rarefied.Rdata")
