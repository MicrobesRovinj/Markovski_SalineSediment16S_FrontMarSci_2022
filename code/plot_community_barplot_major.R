#################################################################################################################
# plot_community_barplot_major.R
#
# A script to plot the relative contribution of the most abundanT taxonomic groups within Desulfobacterota, Gammaproteobacteria, Bacteroidota, Chloroflexi, Planctomycetota, and Campylobacterota.
# Dependencies: results/numerical/community.Rdata
#               data/raw/metadata.csv
#               data/raw/group_colors.csv
# Produces: results/figures/community_barplot_major_1.jpg
#           results/figures/community_barplot_major_2.jpg
#
#################################################################################################################

# Loading taxonomy data
load("results/numerical/community.Rdata")

# Loading colors for each group on the plot
color <- read_tsv("data/raw/group_colors.csv") %>%
  select(- Taxlevel) %>%
  deframe()

# Loading metadata
metadata <- read_tsv("data/raw/metadata.csv")

# Defining station labels for faceting
station_labs <- c("Vegetated", "Nonvegetated")
names(station_labs) <- c("SCy", "SN")

#################################################################################################################
# Generating input data for taxonomy bar plots in a for loop
#################################################################################################################

# Defining taxonomic groups, their corresponding taxonomic levels and threshold values
groups <- tribble(~ taxa, ~ taxlevel, ~ threshold_value,
                  "Desulfobacterota", 2, 2,
                  "Gammaproteobacteria", 3, 2,
                  "Bacteroidota", 2, 2,
                  "Chloroflexi", 2, 1,
                  "Planctomycetota", 2, 1,
                  "Campylobacterota", 2, 1)

for (i in groups$taxa) {
  
  # Selecting taxonomic groups above threshold value
  select <- filter(community,
                   taxlevel == 5 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, paste0("^", i, "$")))$rankID, "."))) %>%
            filter(if_any(6 : ncol(.), ~ . >= filter(groups, taxa == i)$threshold_value))
    
  # Selecting taxonomic groups for plotting
  plot <- filter(community,
                 taxlevel == 5 & str_detect(rankID, paste0("^", filter(community, str_detect(taxon, paste0("^", i, "$")))$rankID, "."))) %>%
    mutate(across(5 : ncol(.), ~ . / sum(.) * 100)) %>%
    filter(rankID %in% select$rankID) %>%
    bind_rows(summarise(., across(everything(.), ~ ifelse(is.numeric(.), 100 - sum(.), paste0("Other_", i))))) %>%
    # Removing last digit from "uncultured" rankID to be used in the next step
    mutate(rankID = if_else(taxon == "uncultured", str_replace(rankID, "\\.\\d+$", ""), rankID)) %>%
    # Removing last two digits from "uncultured_fa" rankID to be used in the next step
    mutate(rankID = if_else(taxon == "uncultured_fa", str_replace(rankID, "(\\.\\d+){2}$", ""), rankID)) %>%
    # Removing last digit from "Unknown_Family" rankID to be used in the next step
    mutate(rankID = if_else(taxon == "Unknown_Family", str_replace(rankID, "\\.\\d+$", ""), rankID))
  
  # Adding information to "uncultured" taxa describing the higher taxonomic levels to which they belong
  uncultured <- select(filter(community, rankID %in% filter(plot, taxon == "uncultured" |
                                                                  taxon == "uncultured_fa" |
                                                                  taxon == "Unknown_Family")$rankID), rankID, taxon) %>%
    rename(rankID_uncultured = rankID, taxon_uncultured = taxon)
  plot <- left_join(plot, uncultured, by = c("rankID" = "rankID_uncultured")) %>%
    mutate(taxon = if_else(taxon == "uncultured", paste0(taxon, "_", taxon_uncultured), taxon)) %>%
    mutate(taxon = if_else(taxon == "uncultured_fa", paste0("uncultured", "_", taxon_uncultured), taxon)) %>%
    mutate(taxon = if_else(taxon == "Unknown_Family", taxon_uncultured, taxon)) %>%
    select(-taxon_uncultured)
  
  # Generating italic names for taxonomic groups in legend
  names <- case_when(i == "Desulfobacterota" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                      str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                                 plot$taxon == "Other_Desulfobacterota" ~ "plain('Other')~italic('Desulfobacterota')",
                                                                      TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     i == "Gammaproteobacteria" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                         str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                                    plot$taxon == "B2M28_fa" ~ "plain('B2M28')",
                                                                                    plot$taxon == "Gammaproteobacteria_Incertae_Sedis" ~ "italic('Gammaproteobacteria Incertae Sedis')",
                                                                                    plot$taxon == "Other_Gammaproteobacteria" ~ "plain('Other')~italic('Gammaproteobacteria')",
                                                                         TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     i == "Bacteroidota" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                  str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                             plot$taxon == "Bacteroidetes_BD2-2" ~ "plain('Bacteroidetes BD2-2')",
                                                                             plot$taxon == "Other_Bacteroidota" ~ "plain('Other')~italic('Bacteroidota')",
                                                                  TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     i == "Chloroflexi" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                 str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                            plot$taxon == "SBR1031_fa" ~ "plain('SBR1031')",
                                                                            plot$taxon == "AB-539-J10" ~ "plain('AB-539-J10')",
                                                                            plot$taxon == "Other_Chloroflexi" ~ "plain('Other')~italic('Chloroflexi')",
                                                                 TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     i == "Planctomycetota" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                     str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                                plot$taxon == "4572-13" ~ "plain('4572-13')",
                                                                                plot$taxon == "SG8-4" ~ "plain('SG8-4')",
                                                                                plot$taxon == "Other_Planctomycetota" ~ "plain('Other')~italic('Planctomycetota')",
                                                                     TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     i == "Campylobacterota" ~ parse(text = case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                                                                      str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                                                                                 plot$taxon == "Other_Campylobacterota" ~ "plain('Other')~italic('Campylobacterota')",
                                                                      TRUE ~ paste0("italic('", plot$taxon, "')"))),
                     TRUE ~ parse(text = paste0("italic('", plot$taxon, "')")))

  # Tidying sequence abundance data
  plot <- plot %>%
    gather(key = "Group", value = "abundance", 6 : ncol(.))
  
  # Joining sequence abundance data and metadata
  plot <- inner_join(metadata, plot, by = c("ID" = "Group")) %>%
    mutate(taxon = factor(taxon, levels = unique(plot$taxon)))

  # Selecting relative abundance of target group in the whole community, tidying the obtained data and joining with metadata
  whole <- community %>%
    filter(taxlevel == filter(groups, taxa == i)$taxlevel) %>%
    gather(key = "Group", value = "abundance", 6 : ncol(.)) %>%
    filter(taxon == i)
  whole <- inner_join(metadata, whole, by = c("ID" = "Group"))

  # Generating mean abundance for different layers
  plot <- plot %>%
    group_by(layer, taxon, station) %>%
    summarise(mean = mean(abundance))

  # Generating mean abundance in the whole community for different layers
  whole <- whole %>%
    group_by(layer, taxon, station) %>%
    summarise(mean = mean(abundance))

  # Renaming and ordering stations and layers, creating a variable for faceting
  plot <- plot %>%
    mutate(station = factor(station, levels = c("SCy", "SN"))) %>%
    mutate(layer = str_replace(layer, "^top$", "Top")) %>%
    mutate(layer = str_replace(layer, "^upper middle$", "Upper \nMiddle")) %>%
    mutate(layer = str_replace(layer, "^lower middle$", "Lower \nMiddle")) %>%
    mutate(layer = str_replace(layer, "^bottom$", "Bottom")) %>%
    mutate(layer = factor(layer, levels = c("Bottom", "Lower \nMiddle", "Upper \nMiddle", "Top"))) %>%
    mutate(facet = i)

  whole <- whole %>%
    mutate(station = factor(station, levels = c("SCy", "SN"))) %>%
    mutate(layer = str_replace(layer, "^top$", "Top")) %>%
    mutate(layer = str_replace(layer, "^upper middle$", "Upper \nMiddle")) %>%
    mutate(layer = str_replace(layer, "^lower middle$", "Lower \nMiddle")) %>%
    mutate(layer = str_replace(layer, "^bottom$", "Bottom")) %>%
    mutate(layer = factor(layer, levels = c("Bottom", "Lower \nMiddle", "Upper \nMiddle", "Top"))) %>%
    mutate(facet = i)

  # Setting names for objects created in loop
  assign(paste0("plot", sep = "_", i), plot)
  assign(paste0("whole", sep = "_", i), whole)
  assign(paste0("names", sep = "_", i), names)
  
  }

#################################################################################################################
# Generating taxonomy bar plots
#################################################################################################################

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
               axis.text = element_text(size = 8, color = "black"), 
               panel.spacing = unit(0.5, "cm"),
               panel.border = element_blank(),
               strip.text.x = element_text(size = 20, face = "bold", hjust = 0.28),
               strip.text.y = element_text(size = 18, face = "bold.italic", margin = margin(l = 0, r = 25)),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background = element_blank(),
               strip.placement = "outside")

# Generating legends in a for loop
for (i in groups$taxa) {

  # Generating the legend
  legend <- get(paste0("plot_", i)) %>%
    mutate(taxon = factor(taxon, levels = rev(levels(taxon)))) %>%
    ggplot() +
    geom_bar(aes(x = layer, y = mean, fill = taxon), stat = "identity", colour = "black", size = 0.3) +
    scale_fill_manual(breaks = levels(get(paste0("plot_", i))$taxon), values = color, labels = get(paste0("names_", i))) +
    theme

  # Extracting legend
  legend <- cowplot::get_legend(legend)

  # Setting names for objects created in loop
  assign(paste0("legend", sep = "_", i), legend)
  
  }

# Generating plot #1
# Preparing data
plot <- full_join(plot_Desulfobacterota, plot_Gammaproteobacteria) %>%
  full_join(., plot_Bacteroidota) %>%
  mutate(taxon = factor(taxon, levels = rev(levels(taxon))))
whole <- full_join(whole_Desulfobacterota, whole_Gammaproteobacteria) %>%
  full_join(., whole_Bacteroidota)
names <- c(names_Desulfobacterota, names_Gammaproteobacteria, names_Bacteroidota)

# Generating plot
p <- plot %>%
  ggplot() +
  geom_bar(aes(x = layer, y = mean, fill = taxon), stat = "identity", colour = "black", size = 0.3) +
  scale_fill_manual(breaks = levels(plot$taxon), values = color, labels = rev(names)) +
  geom_text(data = whole, aes(x = layer, y = 105, label = paste(format(round(mean, 1), nsmall = 1), "%")),
            family = "Times", fontface = "bold", hjust = 0, size = 4) +
  labs(x = "Layer", y = "%") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.215)), breaks = seq(0, 100, by = 10)) +
  facet_grid(factor(facet, levels = c("Desulfobacterota","Gammaproteobacteria","Bacteroidota")) ~ station, labeller = labeller(station = station_labs), switch = "y") +
  coord_capped_flip(bottom = "right") +
  theme

# Removing the legend
p_1 <- p +
  theme(legend.position = "none")

# Generating plot #2
# Preparing data
plot <- full_join(plot_Chloroflexi, plot_Planctomycetota) %>%
  full_join(., plot_Campylobacterota) %>%
  mutate(taxon = factor(taxon, levels = rev(levels(taxon))))
whole <- full_join(whole_Chloroflexi, whole_Planctomycetota) %>%
  full_join(., whole_Campylobacterota)
names <- c(names_Chloroflexi, names_Planctomycetota, names_Campylobacterota)

# Generating plot
p <- plot %>%
  ggplot() +
  geom_bar(aes(x = layer, y = mean, fill = taxon), stat = "identity", colour = "black", size = 0.3) +
  scale_fill_manual(breaks = levels(plot$taxon), values = color, labels = rev(names)) +
  geom_text(data = whole, aes(x = layer, y = 105, label = paste(format(round(mean, 1), nsmall = 1), "%")),
            family = "Times", fontface = "bold", hjust = 0, size = 4) +
  labs(x = "Layer", y = "%") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.215)), breaks = seq(0, 100, by = 10)) +
  facet_grid(factor(facet, levels = c("Chloroflexi","Planctomycetota","Campylobacterota")) ~ station, labeller = labeller(station = station_labs), switch = "y") +
  coord_capped_flip(bottom = "right") +
  theme

# Removing legend
p_2 <- p +
  theme(legend.position = "none")

# Combining generated plots and legends
cowplot_1 <- cowplot::ggdraw() +
  cowplot::draw_plot(p_1, width = 0.755) +
  cowplot::draw_plot(legend_Desulfobacterota, x = 0.745, y = - 0.1221) +
  cowplot::draw_plot(legend_Gammaproteobacteria, x = 0.745, y = - 0.4015) +
  cowplot::draw_plot(legend_Bacteroidota, x = 0.745, y = - 0.7525)

cowplot_2 <- cowplot::ggdraw() +
  cowplot::draw_plot(p_2, width = 0.755) +
  cowplot::draw_plot(legend_Chloroflexi, x = 0.745, y = - 0.121325) +
  cowplot::draw_plot(legend_Planctomycetota, x = 0.745, y = - 0.4975) +
  cowplot::draw_plot(legend_Campylobacterota, x = 0.745, y = - 0.84785)

# Saving generated plots
ggsave("results/figures/community_barplot_major_1.jpg", cowplot_1, width = 297, height = 210, units="mm", bg = "white")
ggsave("results/figures/community_barplot_major_2.jpg", cowplot_2, width = 297, height = 210, units = "mm", bg = "white")
