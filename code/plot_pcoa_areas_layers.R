#################################################################################################################
# plot_pcoa_areas_layers.R
# 
# A script to generate PCoA figures for different areas and layers.
# Dependencies: results/numerical/rarefied.Rdata
#               data/raw/metadata.csv
# Produces: results/figures/pcoa_figure_areas_layers.jpg
#
#################################################################################################################

# Loading rarefied community data
load(file = "results/numerical/rarefied.Rdata")

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv") %>%
  mutate(decay_roots = if_else(decay_roots == "after", "Decay of Roots\nand Rhizomes", decay_roots))

# Joining metadata and rarefied community data
rarefied_metadata <- inner_join(metadata, rarefied, by = c("ID" = "Group"))

# Setting function for number of decimal places
scaleFUN <- function(x) sprintf("%.2f", x)

# Generating a common theme for plots
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               axis.line = element_blank(),
               axis.text = element_text(size = 8, color = "black"),
               axis.title = element_text(size = 10, color = "black"),
               plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
               legend.text = element_text(size = 10),
               legend.text.align = 0,
               legend.key = element_rect(fill = NA),
               legend.key.width = unit(0, "cm"),
               legend.key.height = unit(0.60, "cm"),
               legend.margin = margin(0, 0, 0, 0),
               legend.position = "right",
               legend.justification = c("left", "bottom"))

# Defining colours for each layer
colour_layers <- c("Top" = "#33A02C",
                   "Upper Middle" = "#B2DF8A",
                   "Lower Middle" = "#FFB90F",
                   "Bottom" = "#8B6508")

# Defining colours for each month
colours_months <- c("February" = "#A6CEE3",
                    "March" = "#1F78B4",
                    "April" = "#B2DF8A",
                    "May" = "#33A02C",
                    "June" = "#FB9A99",
                    "July" = "#E31A1C" ,
                    "August" = "#FDBF6F",
                    "September" = "#FF7F00",
                    "October" =  "#CAB2D6",
                    "November" = "#6A3D9A",
                    "December" = "#B15928")

# Defining shapes for each year
shapes_years <- c("2017" = 23,
                  "2018" = 24,
                  "before" = NULL,
                  "Decay of Roots\nand Rhizomes" = 3)

#################################################################################################################
# Generating the PCoA plot for each area and layer in for loops
#################################################################################################################

# Generating the PCoA plot for each area in a for loop
areas <- c("SCy", "SN")

for (i in areas) {

  # Selecting samples for plotting
  rarefied_metadata_select <- rarefied_metadata %>%
    filter(station == i)
  
  # Generating PCoA data
  spe_bray <- rarefied_metadata_select %>%
    column_to_rownames("ID") %>%
    select(starts_with("Otu")) %>%
    select(where(~ sum(.) !=  0)) %>%
    vegdist(method = "bray")
  
  spe_b_pcoa <- wcmdscale(spe_bray, k = (nrow(rarefied_metadata_select) - 1), eig = TRUE)
  
  coordinates <- spe_b_pcoa %>%
    scores(choices = c(1, 2)) %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    rename(A1 = Dim1, A2 = Dim2) %>%
    as_tibble()

  # Joining generated PCoA data and metadata, renaming layers for plotting
  coordinates <- inner_join(metadata, coordinates, by = c("ID" = "ID")) %>%
    mutate(layer = str_replace(layer, "^top$", "Top")) %>%
    mutate(layer = str_replace(layer, "^upper middle$", "Upper Middle")) %>%
    mutate(layer = str_replace(layer, "^lower middle$", "Lower Middle")) %>%
    mutate(layer = str_replace(layer, "^bottom$", "Bottom"))
  
  # Generating the PCoA plot
  p <- ggplot() +
    geom_point(data = coordinates, aes(x = A1, y = A2, fill = layer), shape = 21, size = 2, stroke = 0.5) +
    scale_fill_manual(name = NULL,
                      breaks = names(colour_layers),
                      values = colour_layers) +
    labs(x = paste0("PCoA I (", format(round(spe_b_pcoa$eig[1] / sum(spe_b_pcoa$eig) * 100, digits = 2), nsmall = 2), " %)"), 
         y = paste0("PCoA II (", format(round(spe_b_pcoa$eig[2] / sum(spe_b_pcoa$eig) * 100, digits = 2), nsmall = 2), " %)")) +
    scale_x_continuous(labels = scaleFUN) +
    scale_y_continuous(labels = scaleFUN, breaks = seq(-0.40, 0.40, by = 0.20)) +
    coord_fixed(ratio = 1, xlim = c(-0.50, 0.50), ylim = c(-0.35, 0.45)) +
    theme +
    guides(fill = guide_legend(override.aes = list(size = 3.5)))
  
  # Extracting legend
  legend_1 <- cowplot::get_legend(p)
  
  # Removing legend from plot
  p <- p +
    theme(legend.position = "none")
  
  # Setting name for each plot
  assign(paste0("p", sep = "_", i), p)
  
  # Generating the PCoA plot for each layer and area in a for loop
  layers <- c("top", "upper middle", "lower middle", "bottom")
  
  for (j in layers) {
    
    # Selecting samples for plotting
    rarefied_metadata_select <- rarefied_metadata %>%
      filter(station == i) %>%
      filter(layer == j)
    
    # Generating PCoA data
    spe_bray <- rarefied_metadata_select %>%
      column_to_rownames("ID") %>%
      select(starts_with("Otu")) %>%
      select(where(~ sum(.) != 0)) %>%
      vegdist(method = "bray")
   
     spe_b_pcoa <- wcmdscale(spe_bray, k = (nrow(rarefied_metadata_select) - 1), eig = TRUE)
    
     coordinates <- spe_b_pcoa %>%
       scores(choices = c(1, 2)) %>%
       as.data.frame() %>%
       rownames_to_column("ID") %>%
       rename(A1 = Dim1, A2 = Dim2) %>%
       as_tibble()
    
     # Joining generated PCoA data and metadata, customising dates and years for plotting
     coordinates <- inner_join(metadata, coordinates, by = c("ID" = "ID")) %>%
      mutate(date = as.Date(date, "%d.%m.%Y")) %>%
      mutate(date = format(date, "%B %Y")) %>%
      separate(date, c("month", "year"), sep = " ")
    
    # Generating the PCoA plot
    p <- ggplot() +
      geom_point(data = coordinates, aes(x = A1, y = A2, fill = month, shape = year), size = 2, stroke = 0.5) +
      geom_point(data = coordinates, aes(x = A1, y = A2, shape = decay_roots), size = 3, stroke = 0.3) +
      scale_fill_manual(name = NULL,
                        breaks = names(colours_months),
                        values = colours_months) +
      scale_shape_manual(name = NULL,
                         breaks = names(shapes_years),
                         values = shapes_years) +
      labs(x = paste0("PCoA I (", format(round(spe_b_pcoa$eig[1] / sum(spe_b_pcoa$eig) * 100, digits = 2), nsmall = 2), " %)"), 
           y = paste0("PCoA II (", format(round(spe_b_pcoa$eig[2] / sum(spe_b_pcoa$eig) * 100, digits = 2), nsmall = 2), " %)")) +
      scale_x_continuous(labels = scaleFUN) +
      scale_y_continuous(labels = scaleFUN, breaks = seq(-0.40, 0.40, by = 0.20)) +
      coord_fixed(ratio = 1, xlim = c(-0.50, 0.50), ylim = c(-0.35, 0.45)) +
      theme +
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 3.5)), shape = guide_legend(override.aes = list(size = 3.5)))
    
    # Extracting legend
    legend_2 <- cowplot::get_legend(p)
    
    # Removing legend from plot
    p <- p +
      theme(legend.position = "none")
    
    # Setting name for each plot
    assign(paste0("p", sep = "_", i, sep = "_", j), p)
    
  }
}

# Combining plots and saving
p <- cowplot::plot_grid(p_SCy, p_SN,
                        p_SCy_top, p_SN_top,
                        `p_SCy_upper middle`, `p_SN_upper middle`,
                        `p_SCy_lower middle`, `p_SN_lower middle`,
                        p_SCy_bottom, p_SN_bottom,
                        nrow = 5, ncol = 2)

p <- cowplot::ggdraw() +
  theme(plot.background = element_rect(fill = "#ffffff", color = NA)) +
  cowplot::draw_label("Vegetated", x = 0.280, y = 0.975, hjust = 0.5,  fontfamily = "Times", fontface = "bold", size = 16) +
  cowplot::draw_label("Nonvegetated", x = 0.675, y = 0.975, hjust = 0.5,  fontfamily = "Times", fontface = "bold", size = 16) +
  cowplot::draw_label("All Layers", x = 0.035, y = 0.869, vjust = 0.5, hjust = 0.5, angle = 90, fontfamily = "Times", fontface = "bold", size = 16) +
  cowplot::draw_label("Top", x = 0.035, y = 0.679, vjust = 0.5, hjust = 0.5, angle = 90, fontfamily = "Times", fontface = "bold", size = 16) +
  cowplot::draw_label("Upper Middle", x = 0.035, y = 0.489, vjust = 0.5, hjust = 0.5, angle = 90, fontfamily = "Times", fontface = "bold", size = 16) +
  cowplot::draw_label("Lower Middle", x = 0.035, y = 0.299, vjust = 0.5, hjust = 0.5, angle = 90, fontfamily = "Times", fontface = "bold", size = 16) +
  cowplot::draw_label("Bottom", x = 0.035, y = 0.109, vjust = 0.5, hjust = 0.5, angle = 90, fontfamily = "Times", fontface = "bold", size = 16) +
  cowplot::draw_plot(p, x = 0.05, y = 0, width = 0.80, height = 0.95) +
  cowplot::draw_plot(legend_1, x = 0.84, y = 0.80, width = 0, height = 1) +
  cowplot::draw_plot(legend_2, x = 0.84, y = 0.04, width = 0, height = 1)

ggsave("results/figures/pcoa_figure_areas_layers.jpg", p, width = 195, height = 297, units = "mm")
