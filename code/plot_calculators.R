#################################################################################################################
# plot_calculators.R
#
# A script to plot richness estimators and diversity indices by station and layer.
# Dependencies: results/numerical/estimators_indices_metadata.Rdata
# Produces: results/figures/calculators.jpg
#
#################################################################################################################

# Loading calculated estimators and indices
load("results/numerical/estimators_indices_metadata.Rdata")

# Tidying data for plotting
estimators_indices_metadata <- estimators_indices_metadata %>%
  gather("S.obs", "S.chao1", "S.ACE", "eshannon", "invsimpson", key = "estimator_index", value = "value") 

# Formatting layer, estimator and index names
estimators_indices_metadata <- estimators_indices_metadata %>%
  mutate(layer = str_replace(layer, "^top$", "Top")) %>%
  mutate(layer = str_replace(layer, "^upper middle$", "Upper \nMiddle")) %>%
  mutate(layer = str_replace(layer, "^lower middle$", "Lower \nMiddle")) %>%
  mutate(layer = str_replace(layer, "^bottom$", "Bottom")) %>%
  mutate(layer = factor(layer, levels = c("Top", "Upper \nMiddle", "Lower \nMiddle", "Bottom"))) %>%
  mutate(estimator_index = str_replace(estimator_index, "^S.obs", "Observed Number \nof OTUs")) %>%
  mutate(estimator_index = str_replace(estimator_index, "^S.chao1", "Chao1")) %>%
  mutate(estimator_index = str_replace(estimator_index, "^S.ACE", "ACE")) %>%
  mutate(estimator_index = str_replace(estimator_index, "^eshannon", "Exponential \nShannon")) %>%
  mutate(estimator_index = str_replace(estimator_index, "^invsimpson", "Inverse \nSimpson")) %>%
  mutate(estimator_index = factor(estimator_index, levels = c("Observed Number \nof OTUs", "Chao1", "ACE", "Exponential \nShannon", "Inverse \nSimpson")))

# Generating a common theme for plots
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               axis.line = element_line(colour = "black", lineend = "square", size = 0.2),
               axis.text = element_text(size = 8, color = "black"), 
               axis.title = element_text(size = 10, color = "black"),
               axis.text.x = element_text(size = 10, angle = 0, vjust = 1, hjust = 0.5),
               axis.text.y = element_text(size = 10, angle = 0),
               axis.title.x = element_text(size = 12),
               axis.title.y = element_text(vjust = -18, size = 12),
               plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"),
               legend.position = "none",
               panel.spacing = unit(1, "cm"),
               panel.border = element_blank(),
               strip.text.x = element_text(size = 16, face = "bold", vjust = 2),
               strip.text.y = element_text(size = 16, face = "bold", margin = margin(l = 0, r = 20)),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background = element_blank(),
               strip.placement = "outside")

# Renaming station labels
station_labs <- c("Vegetated", "Nonvegetated")
names(station_labs) <- c("SCy", "SN")

# Generating a function to set custom axis breaks
custom_breaks <- function(x) {
  unlist(case_when(max(x) <=   600                   ~ list(seq(0, 600,  by = 100)),
                   max(x) >   1000 & max(x) <=  1800 ~ list(seq(300, 1800, by = 300)),
                   max(x) >   3000 & max(x) <=  4500 ~ list(seq(1500, 4500, by = 500)),
                   max(x) >  10000                   ~ list(seq(2000, 14000, by = 2000)),
                   TRUE                              ~ list(seq(0, 0, by = 0))))
}

# Generating a function to set custom axis limits
custom_limits <- function(x) {
  case_when(max(x) <=  600                  ~ c(0, 600),
            max(x) >  1000 & max(x) <= 1800 ~ c(300, 1800),
            max(x) >  3000 & max(x) <= 4500 ~ c(1500, 4500),
            max(x) > 10000                  ~ c(2000, 14000),
            TRUE                            ~ c(0, 0))
}

# Adding statistic annotations above box plots
statistics_annotation <- tibble(layer = rep(c("Top", "Upper \nMiddle", "Lower \nMiddle", "Bottom"), 10),
                                value = rep(c(rep(4250, 4),
                                              rep(13000, 4),
                                              rep(13000, 4),
                                              rep(1650, 4),
                                              rep(550, 4)), 2),
                                station = c(rep("SCy", 20), rep("SN", 20)), 
                                estimator_index = rep(c(rep("Observed Number \nof OTUs", 4),
                                                        rep("Chao1", 4),
                                                        rep("ACE", 4),
                                                        rep("Exponential \nShannon", 4),
                                                        rep("Inverse \nSimpson", 4)), 2),
                                label = c(c("a", "ab", "ab", "b"),
                                          rep("a", 8),
                                          c("a", "ab", "bc", "c"),
                                          c("a", "a", "ab", "b"),
                                          rep("a", 4 * 5))) %>%
  mutate(estimator_index = factor(estimator_index, levels = c("Observed Number \nof OTUs", "Chao1", "ACE", 
                                                              "Exponential \nShannon", "Inverse \nSimpson")))

# Generating the plot
p <- estimators_indices_metadata %>%
  ggplot(aes(x = layer, y = value)) +
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(0.5)) +
  geom_boxplot(width = 0.5, fill = "gray90") +
  geom_text(data = statistics_annotation, aes(x = layer, y = value, label = label), family = "Times", fontface = "bold") +
  scale_y_continuous(breaks = custom_breaks, limits = custom_limits, expand = c(0, 0)) +
  labs(x = "Layer", y = "Number of OTUs") +
  facet_rep_grid(estimator_index ~ station, labeller = labeller(station = station_labs), switch = "y", scales = "free_y") +
  theme

# Saving the generated plot
ggsave("results/figures/calculators.jpg", p, width = 180, height = 297, units = "mm")
