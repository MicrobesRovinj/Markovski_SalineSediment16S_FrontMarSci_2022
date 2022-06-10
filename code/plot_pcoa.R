#################################################################################################################
# plot_pcoa.R
# 
# A script to generate the PCoA figure.
# Dependencies: results/numerical/rarefied.Rdata
#               data/raw/metadata.csv
# Produces: results/figures/pcoa_figure.jpg
#
#################################################################################################################

# Loading rarefied community data
load(file = "results/numerical/rarefied.Rdata")

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv")

# Joining metadata and rarefied community data
rarefied_metadata <- inner_join(metadata, rarefied, by = c("ID" = "Group"))

# Generating a common theme for plots
theme <- theme(text = element_text(family = "Times"),
               line = element_line(color = "black"),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               axis.line = element_blank(),
               axis.text = element_text(size = 12, color = "black"),
               axis.title = element_text(size = 14, color = "black"),
               plot.margin = unit(c(5.5, 16.5, 5.5, 16.5), "pt"),
               legend.text = element_text(size = 14, margin = margin(r = 0.2, unit = "cm")),
               legend.text.align = 0,
               legend.key = element_rect(fill = NA),
               legend.key.width = unit(0, "cm"),
               legend.key.height = unit(0.60, "cm"),
               legend.margin = margin(-5, 0, 0, 0),
               legend.spacing.y = unit(0.3, "cm"),
               legend.position = c(0.02, 0.02),
               legend.justification = c("left", "bottom"))

# Generating PCoA data
spe_bray <- rarefied_metadata %>%
  column_to_rownames("ID") %>%
  select(starts_with("Otu")) %>%
  select(where(~ sum(.) !=  0)) %>%
  vegdist(method = "bray")

spe_b_pcoa <- wcmdscale(spe_bray, k = (nrow(rarefied_metadata) - 1), eig = TRUE)

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

# Setting function for number of decimal places
scaleFUN <- function(x) sprintf("%.2f", x)

# Generating the PCoA plot
p <- ggplot() +
  geom_point(data = coordinates, aes(x = A1, y = A2, fill = layer, shape = station), size = 5, stroke = 0.5) +
  scale_fill_manual(name = NULL,
                    breaks = c("Top",
                               "Upper Middle",
                               "Lower Middle",
                               "Bottom"),
                    values = c("Top" = "#33A02C",
                               "Upper Middle" = "#B2DF8A",
                               "Lower Middle" = "#FFB90F",
                               "Bottom" = "#8B6508")) +
  scale_shape_manual(name = NULL,
                     values = c(21, 24),
                     breaks = c("SCy", "SN"),
                     labels = c(parse(text = "Vegetated"), "Nonvegetated")) +
  labs(x = paste0("PCoA I (", format(round(spe_b_pcoa$eig[1] / sum(spe_b_pcoa$eig) * 100, digits = 2), nsmall = 2), " %)"),
       y = paste0("PCoA II (", format(round(spe_b_pcoa$eig[2] / sum(spe_b_pcoa$eig) * 100, digits = 2), nsmall = 2), " %)")) +
  scale_x_continuous(labels = scaleFUN, breaks = seq(-0.4, 0.4, by = 0.2)) +
  scale_y_continuous(labels = scaleFUN, breaks = seq(-0.4, 0.4, by = 0.2)) +
  coord_fixed(ratio = 1, xlim = c(-0.45, 0.45), ylim = c(-0.4, 0.4)) +
  theme +
  guides(fill = guide_legend(override.aes = list(shape = 22)))

# Saving
ggsave("results/figures/pcoa_figure.jpg", p, width = 297, height = 210, units = "mm")
