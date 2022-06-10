#################################################################################################################
# plot_community_barplot_month.R
#
# A script to plot the community structure of each sample showing community temporal dynamics.
# Dependencies: results/numerical/community.Rdata
#               data/raw/metadata.csv
#               data/raw/group_colors.csv
# Produces: results/figures/community_barplot_month.jpg
#
#################################################################################################################

# Loading taxonomy data
load("results/numerical/community.Rdata")

# Selecting taxonomic groups for plotting
plot <- filter(community,
               taxlevel == 2 |
               (taxlevel == 3 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, "^Proteobacteria$"))$rankID, ".")))) %>%
  filter(if_any(6 : ncol(.), ~ . >= 3)) %>%
  mutate(across(5 : ncol(.), ~ case_when(taxon == "Proteobacteria" ~ . - sum(.[taxlevel == 3 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, "^Proteobacteria$"))$rankID, "."))]), TRUE ~ .))) %>%
  mutate(taxon = str_replace(taxon, "^Proteobacteria$", "Other_Proteobacteria")) %>%
  mutate(taxon = str_replace(taxon, "unknown_unclassified", "No_Relative")) %>%
  filter(if_any(6 : ncol(.), ~ . >= 3)) %>%
  bind_rows(summarise(., across(everything(.), ~ ifelse(is.numeric(.), 100 - sum(.), paste("Other"))))) %>%
  arrange(taxon %in% "No_Relative")

# Generating italic names for taxonomic groups in the legend
names <- parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                           plot$taxon == "Sva0485" ~ "plain('Sva0485')",
                                           plot$taxon == "Other_Bacteria" ~ "plain('Other')~italic('Bacteria')",
                                           plot$taxon == "No_Relative" ~ "plain('No Relative')",
                                           plot$taxon == "Other" ~ "plain('Other')",
                                TRUE ~ paste0("italic('", plot$taxon, "')")))

# Tidying sequence abundance data
plot <- plot %>%
  gather(key = "Group", value = "abundance", 6 : ncol(.))

# Loading metadata
metadata <- read_tsv("data/raw/metadata.csv")

# Joining sequence abundance data and metadata
Sys.setlocale(locale = "en_GB.utf8")
plot <- inner_join(metadata, plot, by = c("ID" = "Group")) %>%
  mutate(date = as.Date(date, "%d.%m.%Y")) %>%
  mutate(taxon = factor(taxon, levels = unique(plot$taxon))) %>%
  mutate(label = factor(label, levels = metadata$label)) %>%
  mutate(station = factor(station, levels = c("SCy", "SN")))

# Loading colors for each group on the plot
color <- read_tsv("data/raw/group_colors.csv") %>%
  select(-Taxlevel) %>%
  deframe()

# Defining labels for stations
station_labs <- c("Vegetated", "Nonvegetated")
names(station_labs) <- c("SCy", "SN")

# Renaming and ordering layers
plot <- plot %>%
  mutate(layer = str_replace(layer, "^top$", "Top")) %>%
  mutate(layer = str_replace(layer, "^upper middle$", "Upper Middle")) %>%
  mutate(layer = str_replace(layer, "^lower middle$", "Lower Middle")) %>%
  mutate(layer = str_replace(layer, "^bottom$", "Bottom")) %>%
  mutate(layer = factor(layer, levels = c("Top", "Upper Middle", "Lower Middle", "Bottom")))

# Generating the theme for ggplot
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.background = element_blank(),
               plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
               axis.text.x = element_text(size = 10, angle = 90, vjust = 1.3, hjust = 1),
               axis.text.y = element_text(size = 10, angle = 0),
               axis.title = element_text(size = 12, color = "black"),
               axis.title.x = element_text(hjust = 0.5),
               axis.title.y = element_text(vjust = -15),
               axis.line.x = element_line(),
               axis.line.y = element_line(),
               legend.position = "right",
               legend.text = element_text(size = 10),
               legend.text.align = 0,
               legend.spacing.x = unit(0.2, "cm"),
               legend.justification = c("left", "bottom"),
               legend.key = element_rect(fill = NA),
               legend.key.size = unit(0.45, "cm"),
               legend.box = "vertical",
               legend.title = element_blank(),
               panel.grid = element_blank(),
               axis.text = element_text(size = 10, color = "black"), 
               panel.spacing = unit(0.8, "cm"),
               panel.border = element_blank(),
               strip.text.x = element_text(size = 20, face = "bold", vjust = 2),
               strip.text.y = element_text(size = 20, face = "bold", margin = margin(l = 0, r = 30)),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background = element_blank(),
               strip.placement = "outside")

# Generating the plot
p <- ggplot(plot) +
  geom_bar(aes(x = date, y = abundance, fill = taxon), stat = "identity", colour = "black", size = 0.3) +
  scale_fill_manual(breaks = levels(plot$taxon), values = color, labels = names, guide = guide_legend(reverse = FALSE, ncol = 1)) +
  labs(x = "Date", y = "%") +
  scale_y_continuous(breaks = seq(0, 100, by = 10), expand = c(0, 0)) +
  scale_x_date(breaks = seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels = c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                          "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                          "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits = as.Date(c("2017-07-01", "2018-11-01")),
               expand = c(0, 0)) +
  facet_rep_grid(layer ~ station, labeller = labeller(station = station_labs), switch = "y") +
  coord_capped_cart(left = 'both') +
  theme

# Saving the generated plot
ggsave("results/figures/community_barplot_month.jpg", p, width = 210, height = 297, units = "mm", bg = "white")
