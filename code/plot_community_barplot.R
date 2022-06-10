#################################################################################################################
# plot_community_barplot.R
#
# A script to plot the relative contribution of the most abundant bacterial and archaeal taxonomic groups.
# Dependencies: results/numerical/community.Rdata
#               data/raw/metadata.csv
#               data/raw/group_colors.csv
# Produces: results/figures/community_barplot.jpg
#
#################################################################################################################

# Loading taxonomy data
load("results/numerical/community.Rdata")

#################################################################################################################
# Archaea
#################################################################################################################

# Selecting taxonomic groups above threshold value
select_archaea <- filter(community,
                         taxlevel == 2 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, "^Archaea$"))$rankID, "."))) %>%
  filter(if_any(6 : ncol(.), ~ . >= 3))

# Selecting taxonomic groups for plotting
plot_archaea <- filter(community,
                       taxlevel == 2 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, "^Archaea$"))$rankID, "."))) %>%
  mutate(across(5 : ncol(.), ~ . / sum(.) * 100)) %>%
  filter(rankID %in% select_archaea$rankID) %>%
  bind_rows(summarise(., across(everything(.), ~ ifelse(is.numeric(.), 100 - sum(.), paste("Other_Archaea"))))) %>%
  # Removing last digit from "uncultured" rankID to be used in the next step
  mutate(rankID = if_else(taxon == "uncultured", str_replace(rankID, "\\.\\d+$", ""), rankID))

# Adding information to "uncultured" taxa describing the higher taxonomic level to which they belong
uncultured <- select(filter(community, rankID %in% filter(plot_archaea, taxon == "uncultured")$rankID), rankID, taxon) %>%
  rename(rankID_uncultured = rankID, taxon_uncultured = taxon)
plot_archaea <- left_join(plot_archaea, uncultured, by = c("rankID" = "rankID_uncultured")) %>%
  mutate(taxon = if_else(taxon == "uncultured", paste0(taxon, "_", taxon_uncultured), taxon)) %>%
  select(-taxon_uncultured)

# Generating italic names for taxonomic groups in the legend
names_archaea <- parse(text = case_when(str_detect(plot_archaea$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot_archaea$taxon, "uncultured_"), "')"),
                                        str_detect(plot_archaea$taxon, "unclassified") ~ paste0("italic('", str_remove(plot_archaea$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                   plot_archaea$taxon == "Other_Archaea" ~ "plain('Other')~italic('Archaea')",
                                        TRUE ~ paste0("italic('", plot_archaea$taxon, "')")))

# Tidying sequence abundance data
plot_archaea <- plot_archaea %>%
  gather(key = "Group", value = "abundance", 6 : ncol(.))

# Loading metadata
metadata <- read_tsv("data/raw/metadata.csv")

# Joining sequence abundance data and metadata
plot_archaea <- inner_join(metadata, plot_archaea, by = c("ID" = "Group")) %>%
  mutate(taxon = factor(taxon, levels = unique(plot_archaea$taxon)))

# Selecting relative abundance of target group in the whole community, tidying the obtained data and joining with metadata
whole_archaea <- community %>%
  filter(taxlevel == 1) %>%
  gather(key = "Group", value = "abundance", 6 : ncol(.)) %>%
  filter(taxon == "Archaea")
whole_archaea <- inner_join(metadata, whole_archaea, by = c("ID" = "Group"))

# Generating mean abundance for different layers
plot_archaea <- plot_archaea %>%
  group_by(layer, taxon, station) %>%
  summarise(mean = mean(abundance))

# Generating mean abundance in the whole community for different layers
whole_archaea <- whole_archaea %>%
  group_by(layer, taxon, station) %>%
  summarise(mean = mean(abundance))

# Renaming and ordering stations and layers, creating a variable for faceting
plot_archaea <- plot_archaea %>%
  mutate(layer = str_replace(layer, "^top$", "Top")) %>%
  mutate(layer = str_replace(layer, "^upper middle$", "Upper \nMiddle")) %>%
  mutate(layer = str_replace(layer, "^lower middle$", "Lower \nMiddle")) %>%
  mutate(layer = str_replace(layer, "^bottom$", "Bottom")) %>%
  mutate(layer = factor(layer, levels = c("Bottom", "Lower \nMiddle", "Upper \nMiddle", "Top"))) %>%
  mutate(facet = "Archaea")

whole_archaea <- whole_archaea %>%
  mutate(layer = str_replace(layer, "^top$", "Top")) %>%
  mutate(layer = str_replace(layer, "^upper middle$", "Upper \nMiddle")) %>%
  mutate(layer = str_replace(layer, "^lower middle$", "Lower \nMiddle")) %>%
  mutate(layer = str_replace(layer, "^bottom$", "Bottom")) %>%
  mutate(layer = factor(layer, levels = c("Bottom", "Lower \nMiddle", "Upper \nMiddle", "Top"))) %>%
  mutate(facet = "Archaea")

#################################################################################################################
# Bacteria
#################################################################################################################

# Selecting taxonomic groups above threshold value
select_bacteria <- filter(community,
                          (taxlevel == 2 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, "^Bacteria$"))$rankID, "."))) |
                          (taxlevel == 3 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, "^Proteobacteria$"))$rankID, ".")))) %>%
  filter(if_any(6 : ncol(.), ~ . >= 3)) %>%
  mutate(across(5 : ncol(.), ~ case_when(taxon == "Proteobacteria" ~ . - sum(.[taxlevel == 3 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, "^Proteobacteria$"))$rankID, "."))]), TRUE ~ .))) %>%
  mutate(taxon = str_replace(taxon, "^Proteobacteria$", "Other_Proteobacteria")) %>%
  filter(if_any(6 : ncol(.), ~ . >= 3))

# Selecting taxonomic groups for plotting
plot_bacteria <- filter(community,
                        (taxlevel == 2 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, "^Bacteria$"))$rankID, "."))) |
                        (taxlevel == 3 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, "^Proteobacteria$"))$rankID, ".")))) %>%
  mutate(across(5 : ncol(.), ~ case_when(taxon == "Proteobacteria" ~ . - sum(.[taxlevel == 3 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, "^Proteobacteria$"))$rankID, "."))]), TRUE ~ .))) %>%
  mutate(taxon = str_replace(taxon, "^Proteobacteria$", "Other_Proteobacteria")) %>%
  mutate(across(5 : ncol(.), ~ . / sum(.) * 100)) %>%
  filter(rankID %in% select_bacteria$rankID) %>%
  bind_rows(summarise(., across(everything(.), ~ ifelse(is.numeric(.), 100 - sum(.), paste("Other_Bacteria"))))) %>%
  # Removing last digit from "uncultured" rankID to be used in the next step
  mutate(rankID = if_else(taxon == "uncultured", str_replace(rankID, "\\.\\d+$", ""), rankID))

# Adding information to "uncultured" taxa describing the higher taxonomic level to which they belong
uncultured <- select(filter(community, rankID %in% filter(plot_bacteria, taxon == "uncultured")$rankID), rankID, taxon) %>%
  rename(rankID_uncultured = rankID, taxon_uncultured = taxon)
plot_bacteria <- left_join(plot_bacteria, uncultured, by = c("rankID" = "rankID_uncultured")) %>%
  mutate(taxon = if_else(taxon == "uncultured", paste0(taxon, "_", taxon_uncultured), taxon)) %>%
  select(-taxon_uncultured)

# Generating italic names for taxonomic groups in the legend
names_bacteria <- parse(text = case_when(str_detect(plot_bacteria$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot_bacteria$taxon, "uncultured_"), "')"),
                                         str_detect(plot_bacteria$taxon, "unclassified") ~ paste0("italic('", str_remove(plot_bacteria$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                    plot_bacteria$taxon == "Sva0485" ~ "plain('Sva0485')",
                                                    plot_bacteria$taxon == "Other_Bacteria" ~ "plain('Other')~italic('Bacteria')",
                                         TRUE ~ paste0("italic('", plot_bacteria$taxon, "')")))

# Tidying sequence abundance data
plot_bacteria <- plot_bacteria %>%
  gather(key = "Group", value = "abundance", 6 : ncol(.))

# Joining sequence abundance data and metadata
plot_bacteria <- inner_join(metadata, plot_bacteria, by = c("ID" = "Group")) %>%
  mutate(taxon = factor(taxon, levels = unique(plot_bacteria$taxon)))

# Selecting relative abundance of target group in the whole community, tidying the obtained data and joining with metadata
whole_bacteria <- community %>%
  filter(taxlevel == 1) %>%
  gather(key = "Group", value = "abundance", 6 : ncol(.)) %>%
  filter(taxon == "Bacteria")
whole_bacteria <- inner_join(metadata, whole_bacteria, by = c("ID" = "Group"))

# Generating mean abundance for different layers
plot_bacteria <- plot_bacteria %>%
  group_by(layer, taxon, station) %>%
  summarise(mean = mean(abundance))

# Generating mean abundance in the whole community for different layers
whole_bacteria <- whole_bacteria %>%
  group_by(layer, taxon, station) %>%
  summarise(mean = mean(abundance))

# Renaming and ordering stations and layers, creating a variable for faceting
plot_bacteria <- plot_bacteria %>%
  mutate(taxon == factor(taxon, levels = unique(plot_bacteria$taxon))) %>%
  mutate(layer = str_replace(layer, "^top$", "Top")) %>%
  mutate(layer = str_replace(layer, "^upper middle$", "Upper \nMiddle")) %>%
  mutate(layer = str_replace(layer, "^lower middle$", "Lower \nMiddle")) %>%
  mutate(layer = str_replace(layer, "^bottom$", "Bottom")) %>%
  mutate(layer = factor(layer, levels = c("Bottom", "Lower \nMiddle", "Upper \nMiddle", "Top"))) %>%
  mutate(facet = "Bacteria")

whole_bacteria <- whole_bacteria %>%
  mutate(layer = str_replace(layer, "^top$", "Top")) %>%
  mutate(layer = str_replace(layer, "^upper middle$", "Upper \nMiddle")) %>%
  mutate(layer = str_replace(layer, "^lower middle$", "Lower \nMiddle")) %>%
  mutate(layer = str_replace(layer, "^bottom$", "Bottom")) %>%
  mutate(layer = factor(layer, levels = c("Bottom", "Lower \nMiddle", "Upper \nMiddle", "Top"))) %>%
  mutate(facet = "Bacteria")

#################################################################################################################
# Generating a common plot
#################################################################################################################

# Generating named vector for labels
names <- c(names_archaea, names_bacteria)

# Joining archaeal and bacterial data for faceting
plot <- full_join(plot_archaea, plot_bacteria)
plot <- plot %>%
  mutate(taxon = factor(taxon, levels = rev(levels(plot$taxon))))
whole <- full_join(whole_archaea, whole_bacteria)

# Loading colors for each group on the plot
color <- read_tsv("data/raw/group_colors.csv") %>%
  select(-Taxlevel) %>%
  deframe()

# Defining labels for stations and facets
station_labs <- c("Vegetated", "Nonvegetated")
names(station_labs) <- c("SCy", "SN")

facet_labs <- c("Bacteria", "Archaea")
names(facet_labs) <- c("Bacteria", "Archaea")

# Generating the theme for ggplot
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.background = element_blank(),
               plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
               axis.text.x = element_text(size = 10, angle = 0, vjust = 1, hjust = 0.5),
               axis.text.y = element_text(size = 10, angle = 0),
               axis.title = element_text(size = 12, color = "black"),
               axis.title.x = element_text(hjust = 0.45),
               axis.title.y = element_text(vjust = -15),
               axis.line.x = element_line(),
               legend.position = "right",
               legend.text = element_text(size = 12),
               legend.text.align = 0,
               legend.spacing.x = unit(0.2, "cm"),
               legend.justification = c("left", "top"),
               legend.key = element_rect(fill = NA),
               legend.key.size = unit(0.5, "cm"),
               legend.box = "vertical",
               legend.title = element_blank(),
               panel.grid = element_blank(),
               axis.text = element_text(size = 10, color = "black"),
               panel.spacing = unit(0.5, "cm"),
               panel.border = element_blank(),
               strip.text.x = element_text(size = 20, face = "bold", hjust = 0.3),
               strip.text.y = element_text(size = 20, face = "bold.italic", margin = margin(l = 0, r = 20)),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background = element_blank(),
               strip.placement = "outside")

# Generating the plot
p <- ggplot(plot) +
  geom_bar(aes(x = layer, y = mean, fill = taxon), stat = "identity", colour = "black", size = 0.3) +
  scale_fill_manual(breaks = rev(levels(plot$taxon)), values = color, labels = names) +
  geom_text(data = whole, aes(x = layer, y = 105, label = paste(format(round(mean, 1), nsmall = 1), "%")),
            family = "Times", fontface = "bold", hjust = 0, size = 4) +
  labs(x = "Layer", y = "%") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), breaks = seq(0, 100, by = 10)) +
  facet_grid(facet ~ station, labeller = labeller(station = station_labs, facet = facet_labs), switch = "y") +
  coord_capped_flip(bottom = "right") +
  theme

# Removing legend (to be added later)
p <- p +
  theme(legend.position = "none")

# Generating the legend for Archaea
legend_archaea <- plot %>%
  filter(facet == "Archaea") %>%
  mutate(taxon = factor(taxon, levels = unique(.$taxon)))
legend_archaea <- ggplot(legend_archaea) +
  geom_bar(aes(x = layer, y = mean, fill = taxon), stat = "identity", colour = "black", size = 0.3) +
  scale_fill_manual(breaks = levels(legend_archaea$taxon), values = color, labels = names_archaea) +
  theme
legend_archaea <- cowplot::get_legend(legend_archaea)

# Generating the legend for Bacteria
legend_bacteria <- plot %>%
  filter(facet == "Bacteria") %>%
  mutate(taxon = factor(taxon, levels = unique(.$taxon)))
legend_bacteria <- ggplot(legend_bacteria) +
  geom_bar(aes(x = layer, y = mean, fill = taxon), stat = "identity", colour = "black", size = 0.3) +
  scale_fill_manual(breaks = levels(legend_bacteria$taxon), values = color, labels = names_bacteria) +
  theme
legend_bacteria <- cowplot::get_legend(legend_bacteria)

# Combining the plot with the legends
p <- cowplot::ggdraw() +
  cowplot::draw_plot(p, width = 0.755) +
  cowplot::draw_plot(legend_archaea, x = 0.745, y = - 0.34) + 
  cowplot::draw_plot(legend_bacteria, x = 0.745, y = - 0.533) 

# Saving the generated plot
ggsave("results/figures/community_barplot.jpg", p, width = 297, height = 210, units = "mm", bg = "white")
