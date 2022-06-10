#################################################################################################################
# plot_matrix.R
#
# A script to plot Bray-Curtis similarity coefficients between different stations and layers.
# Dependencies: results/numerical/rarefied.Rdata
#               data/raw/metadata.csv
# Produces: results/figures/matrix.jpg
#
#################################################################################################################

# Loading rarefied data
load(file = "results/numerical/rarefied.Rdata")

# Loading metadata
metadata <- read_tsv("data/raw/metadata.csv")

# Joining metadata with OTU/sample data and summing sequences from each environment
rarefied_metadata <- inner_join(rarefied, metadata, by = c("Group" = "ID")) %>%
  group_by(station, layer) %>%
  summarise(across(starts_with("Otu"), sum), .groups = "drop")

# Sorting stations/layers and copying station/layer labels to row names (input for library vegan)
rarefied_metadata <- rarefied_metadata %>%
  unite("station", station:layer, sep = "_", remove = TRUE) %>%
  mutate(station = factor(station, levels = c("SCy_top","SCy_upper middle", "SCy_lower middle", "SCy_bottom",
                                              "SN_top", "SN_upper middle", "SN_lower middle", "SN_bottom"))) %>%
  arrange(station) %>%
  column_to_rownames("station")

# Calculating Bray-Curtis dissimilarity
bray <- rarefied_metadata %>%
  vegdist(method = "bray", binary = FALSE) %>%
  as.matrix()
bray[upper.tri(bray, diag = TRUE)] <- NA
bray <- bray %>%
  as_tibble(.name_repair = "check_unique", rownames = NA) %>%
  rownames_to_column("V1") %>%
  gather(key = "V2", value = "bray", 2 : ncol(.)) %>%
  filter(!is.na(bray))

# Converting dissimilarity to similarity
similarity <- bray %>%
  mutate(bray = 1 - bray)

# Renaming stations/layers for plotting
similarity <- similarity %>%
  mutate(across(c(V1, V2), ~ str_replace(., "SCy_top", "Top Vegetated Layer"))) %>%
  mutate(across(c(V1, V2), ~ str_replace(., "SCy_upper middle", "Upper Middle Vegetated Layer"))) %>%
  mutate(across(c(V1, V2), ~ str_replace(., "SCy_lower middle", "Lower Middle Vegetated Layer"))) %>%
  mutate(across(c(V1, V2), ~ str_replace(., "SCy_bottom", "Bottom Vegetated Layer"))) %>%
  mutate(across(c(V1, V2), ~ str_replace(., "SN_top", "Top Nonvegetated Layer"))) %>%
  mutate(across(c(V1, V2), ~ str_replace(., "SN_upper middle", "Upper Middle Nonvegetated Layer"))) %>%
  mutate(across(c(V1, V2), ~ str_replace(., "SN_lower middle", "Lower Middle Nonvegetated Layer"))) %>%
  mutate(across(c(V1, V2), ~ str_replace(., "SN_bottom", "Bottom Nonvegetated Layer")))

# Arranging factors for plotting
similarity <- similarity %>%
  mutate(V1 = factor(V1, levels = rev(c("Top Vegetated Layer",
                                        "Upper Middle Vegetated Layer",
                                        "Lower Middle Vegetated Layer",
                                        "Bottom Vegetated Layer",
                                        "Top Nonvegetated Layer",
                                        "Upper Middle Nonvegetated Layer",
                                        "Lower Middle Nonvegetated Layer",
                                        "Bottom Nonvegetated Layer")))) %>%
  mutate(V2 = factor(V2, levels = c("Top Vegetated Layer",
                                    "Upper Middle Vegetated Layer",
                                    "Lower Middle Vegetated Layer",
                                    "Bottom Vegetated Layer",
                                    "Top Nonvegetated Layer",
                                    "Upper Middle Nonvegetated Layer",
                                    "Lower Middle Nonvegetated Layer",
                                    "Bottom Nonvegetated Layer")))

# Setting function for the number of decimal places
scaleFUN <- function(x) sprintf("%.2f", x)

# Generating the theme for ggplot
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.grid = element_blank(),
               axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y.left = element_blank(),
               axis.text.x = element_text(size = 15, color = "black", hjust = 0.5, vjust = 0.5, margin = margin(t = -0.2, unit = "cm"), face = "bold"),
               axis.text.y = element_text(size = 15, color = "black", hjust = 0.5, vjust = 0.5, margin = margin(r = -0.2, unit = "cm"), face = "bold"),
               panel.background = element_blank(),
               plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
               legend.position = "none")

# Generating the plot
p <- ggplot(similarity, aes(x = V2, y = V1, fill = bray)) +
  geom_tile(colour = "black", size = 0.5) +
  geom_text(aes(V2, V1, label = scaleFUN(bray)), color = "black", size = 7, family = "Times", fontface = "bold") +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 15)) +
  scale_fill_gradient2(low = "lightgoldenrod", mid = "lightgoldenrod3", high = "lightgoldenrod4", midpoint = 0.5, limit = c(0, 1)) +
  theme

# Saving the generated plot
ggsave("results/figures/matrix.jpg", p, width = 1.75 * 210, height = 297, units = "mm")
