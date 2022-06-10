#################################################################################################################
# plot_rarefaction.R
# 
# A script to plot the rarefaction curve of each sample.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.groups.rarefaction
#               data/raw/metadata.csv
# Produces: results/figures/rarefaction_a.jpg
#           results/figures/rarefaction_b.jpg
#
#################################################################################################################

# Loading input data and selecting values for plotting
rarefaction <- read_tsv(file = "data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.groups.rarefaction") %>%
  select(-contains("lci-"), -contains("hci-")) %>%
  gather(-numsampled, key = sample, value = sobs) %>%
  mutate(sample = str_replace_all(sample, pattern = "0.03-", replacement = "")) %>%
  drop_na()

# Loading metadata
Sys.setlocale(locale = "en_GB.utf8")
metadata <- read_tsv("data/raw/metadata.csv") %>%
  mutate(date = as.Date(date, "%d.%m.%Y")) %>%
  arrange(date) %>%
  mutate(date_raw = date) %>%
  mutate(date = format(date, "%d %B %Y")) %>%
  mutate(date = str_replace(date, "^0", ""))

# Joining metadata and input data
metadata_rarefaction <- inner_join(metadata, rarefaction, by = c("ID" = "sample")) %>%
  mutate(date = factor(date, levels = unique(date))) %>%
  mutate(depth = str_replace(depth, "-", " – "))

# Defining line colour and type
colour_line_type <- tribble(~station, ~depth, ~colour, ~line_type,
                            "SCy", "0 – 1 cm", "#1F78B4", "solid",
                            "SCy", "1 – 2 cm", "#A6CEE3", "solid",
                            "SCy", "2 – 3 cm", "#33A02C", "solid",
                            "SCy", "3 – 4 cm", "#B2DF8A", "solid",
                            "SCy", "4 – 5 cm", "#E31A1C", "solid",
                            "SCy", "5 – 6 cm", "#FB9A99", "solid",
                            "SCy", "7 – 8 cm", "#FF7F00", "solid",
                            "SN", "0 – 1 cm", "#1F78B4", "dotted",
                            "SN", "1 – 2 cm", "#A6CEE3", "dotted",
                            "SN", "2 – 3 cm", "#33A02C", "dotted",
                            "SN", "3 – 4 cm", "#B2DF8A", "dotted",
                            "SN", "4 – 5 cm", "#E31A1C", "dotted",
                            "SN", "5 – 6 cm", "#FB9A99", "dotted",
                            "SN", "7 – 8 cm", "#FF7F00", "dotted")

# Generating a common theme for plots
theme <- theme(text = element_text(family = "Times"),
               line = element_line(colour = "black"),
               panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(), 
               axis.line = element_line(colour = "gray60"),
               axis.ticks = element_line(colour = "gray60"),
               axis.text = element_text(size = 12, colour = "black"), 
               axis.title = element_text(size = 14, colour = "black"),
               plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"), 
               legend.position = "right",
               legend.title = element_blank(),
               legend.text = element_text(size = 10, margin = margin(r = 0.2, unit = "cm")),
               legend.key.width = unit(1.4, "cm"),
               legend.key.height = unit(0.5, "cm"),
               legend.key = element_rect(fill = "white"), 
               legend.justification = c("bottom"),
               legend.text.align = 0,
               legend.spacing.x = unit(0, "cm"),
               plot.title = element_text(size = 16, hjust = 0.5), 
               strip.background = element_blank(),
               strip.text = element_text(size = 14, face = "bold"),
               legend.margin = margin(t = -10, unit = "pt"))

# Setting function for custom axis breaks
custom_breaks <- function(x) {
  unlist(case_when(max(x) <= 8000                    ~ list(seq(0, 8000,  by=1000)),
                   max(x) >  8000  & max(x) <= 10000 ~ list(seq(0, 10000, by=2000)),
                   max(x) >  10000 & max(x) <= 20000 ~ list(seq(0, 20000, by=5000)),
                   TRUE                              ~ list(seq(0, 60000, by=10000))))
}

# Setting function for custom axis limits
custom_limits <- function(x) {
  case_when(max(x) >= 1000  & max(x) < 10000  ~ c(0, ceiling(max(x) / 1000)  * 1000),
            max(x) >= 10000 & max(x) < 100000 ~ c(0, ceiling(max(x) / 10000) * 10000),
            TRUE                              ~ c(0, max(x)))
}

# Generating plots
p1 <- metadata_rarefaction %>%
  filter(date_raw <= "2018-03-26") %>%
  ggplot() +
  geom_line(aes(x = numsampled, y = sobs, colour = depth, linetype = station), size = 1.0) +
  scale_colour_manual(values = set_names(colour_line_type$colour, colour_line_type$depth)) +
  scale_linetype_manual(values = set_names(colour_line_type$line_type, colour_line_type$station),
                        labels = c("SCy" = "Vegetated", "SN" = "Nonvegetated")) +
  scale_x_continuous(breaks = custom_breaks, limits = custom_limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = custom_breaks, limits = custom_limits, expand = c(0, 0)) +
  labs(x = "Number of Sequences", y = "Number of OTUs") +
  theme +
  facet_wrap(~ date, nrow = 5, ncol = 2, scales = "free")

p2 <- metadata_rarefaction %>%
  filter(date_raw > "2018-03-26") %>%
  ggplot() +
  geom_line(aes(x = numsampled, y = sobs, colour = depth, linetype = station), size = 1.0) +
  scale_colour_manual(values = set_names(colour_line_type$colour, colour_line_type$depth)) +
  scale_linetype_manual(values = set_names(colour_line_type$line_type, colour_line_type$station),
                        labels = c("SCy" = "Vegetated", "SN" = "Nonvegetated")) +
  scale_x_continuous(breaks = custom_breaks, limits = custom_limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = custom_breaks, limits = custom_limits, expand = c(0, 0)) +
  labs(x = "Number of Sequences", y = "Number of OTUs") +
  theme +
  facet_wrap(~ date, nrow = 5, ncol = 2, scales = "free")

# Saving
ggsave("results/figures/rarefaction_a.jpg", p1, width = 210, height = 297, units = "mm")
ggsave("results/figures/rarefaction_b.jpg", p2, width = 210, height = 297*4/5, units = "mm")
