#################################################################################################################
# plot_estimators_per_month.R
#
# A script to plot richness estimators and diversity indices for each sampling point.
# Dependencies: results/numerical/estimators_indices_metadata.Rdata
# Produces: results/figures/diversity_indices_month.jpg
#
#################################################################################################################

# Loading calculated estimators and indices
load("results/numerical/estimators_indices_metadata.Rdata")

# Tidying data for plotting
Sys.setlocale(locale = "en_GB.utf8")
estimators_indices_metadata <- estimators_indices_metadata %>%
  gather("S.obs", "S.chao1", "S.ACE", "eshannon", "invsimpson", key = "estimator_index", value = "value")

# Formatting layer, estimator and index names
estimators_indices_metadata <- estimators_indices_metadata %>%
  mutate(layer = str_replace(layer, "^top$", "Top")) %>%
  mutate(layer = str_replace(layer, "^upper middle$", "Upper Middle")) %>%
  mutate(layer = str_replace(layer, "^lower middle$", "Lower Middle")) %>%
  mutate(layer = str_replace(layer, "^bottom$", "Bottom")) %>%
  mutate(layer = factor(layer, levels = c("Top", "Upper Middle", "Lower Middle", "Bottom"))) %>%
  mutate(estimator_index = str_replace(estimator_index, "^S.obs", "Observed Number \nof OTUs")) %>%
  mutate(estimator_index = str_replace(estimator_index, "^S.chao1", "Chao1")) %>%
  mutate(estimator_index = str_replace(estimator_index, "^S.ACE", "ACE")) %>%
  mutate(estimator_index = str_replace(estimator_index, "^eshannon", "Exponential Shannon")) %>%
  mutate(estimator_index = str_replace(estimator_index, "^invsimpson", "Inverse Simpson")) %>%
  mutate(estimator_index = factor(estimator_index, levels = c("Observed Number \nof OTUs", "Chao1", "ACE", "Exponential Shannon", "Inverse Simpson")))

# Defining line types, dot shapes, and dot fill colors
lines_p1 <- c("Observed Number \nof OTUs" = "dotted", "Chao1" = "solid", "ACE" = "dotted")
lines_p2 <- c("Exponential Shannon" = "solid", "Inverse Simpson" = "dotted")
shapes_p1 <- c("Observed Number \nof OTUs" = 21, "Chao1" = 23, "ACE" = 25)
shapes_p2 <- c("Exponential Shannon" = 21, "Inverse Simpson" = 24)
fills_p1 <- c("Observed Number \nof OTUs" = "white", "Chao1" = "black", "ACE" = "white")
fills_p2 <- c("Exponential Shannon" = "black", "Inverse Simpson" = "white")

# Generating a common theme for plots
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               axis.line = element_line(colour = "black", lineend = "square", size = 0.2),
               axis.text = element_text(size = 8, color = "black"),
               axis.title = element_text(size = 12, color = "black"),
               axis.title.y = element_text(vjust = -13),
               axis.text.x = element_text(size = 10, angle = 90, vjust = 1.25, hjust = 0.95),
               axis.text.y = element_text(size = 10, angle = 0),
               plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
               legend.text = element_text(size = 10, margin = margin(r = 5.5, unit = "pt")), 
               legend.key.width = unit(1, "cm"),
               legend.key.height = unit(0.75, "cm"),
               legend.key = element_rect(fill = "white"),
               legend.text.align = 0, 
               legend.margin = margin(5.5, 5.5, 5.5, 5.5),
               legend.position = "right",
               legend.justification = c("right", "bottom"),
               legend.title = element_blank(),
               plot.title = element_text(size = 14, hjust = 0.5),
               panel.spacing = unit(1, "cm"),
               panel.border = element_blank(),
               strip.text.x = element_text(size = 16, face = "bold", vjust = 2),
               strip.text.y = element_text(size = 16, face = "bold", vjust = 0, margin = margin(r = 25)),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background = element_blank(),
               strip.placement = "outside")

# Defining names for stations and layers to be used for faceting
station_labs <- c("Vegetated", "Nonvegetated")
names(station_labs) <- c("SCy", "SN")
layer_labs <- c("Top", "Upper Middle", "Lower Middle", "Bottom")
names(layer_labs) <- c("Top", "Upper Middle", "Lower Middle", "Bottom")

# Plotting richness estimators
p1 <- estimators_indices_metadata %>%
  filter(estimator_index == "Observed Number \nof OTUs" | estimator_index == "Chao1" | estimator_index == "ACE") %>%
  ggplot(aes(x = date, y = value, linetype = estimator_index, shape = estimator_index,
             fill = estimator_index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines_p1) +
  scale_shape_manual(values = shapes_p1) +
  scale_fill_manual(values = fills_p1) +
  scale_y_continuous(breaks = seq(2000, 12000, by = 2000), limits = c(1000, 12000), expand = c(0, 0)) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                          "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                          "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-06-15", "2018-11-01")),
               expand = c(0, 0)) +
  labs(x = "Date", y = "Number of OTUs") +
  facet_rep_grid(layer ~ station, labeller = labeller(station = station_labs, layer = layer_labs), switch = "y", repeat.tick.labels = FALSE) +
  theme

# Saving generated plot
ggsave("results/figures/estimators_months.jpg", p1, width = 210, height = 297, units = "mm")

# Plotting diversity indices
p2 <- estimators_indices_metadata %>%
  filter(estimator_index == "Exponential Shannon" | estimator_index == "Inverse Simpson") %>%
  ggplot(aes(x = date, y = value, linetype = estimator_index, shape = estimator_index,
             fill = estimator_index)) +
  geom_line() +
  geom_point(size = 3) +
  scale_linetype_manual(values = lines_p2) +
  scale_shape_manual(values = shapes_p2) +
  scale_fill_manual(values = fills_p2) +
  scale_y_continuous(breaks = seq(0, 1500, by = 250), limits = c(0, 1500), expand = c(0, 0)) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                          "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                          "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-06-15", "2018-11-01")),
               expand = c(0, 0)) +
  labs(x = "Date", y = "Number of OTUs") +
  facet_rep_grid(layer ~ station, labeller = labeller(station = station_labs, layer = layer_labs), switch = "y", repeat.tick.labels = FALSE) +
  theme

# Saving generated plot
ggsave("results/figures/diversity_indices_month.jpg", p2, width = 210, height = 297, units = "mm")
