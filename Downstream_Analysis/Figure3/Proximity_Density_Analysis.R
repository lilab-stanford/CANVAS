rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
})

# =========================================================
# 1. Define paths
# =========================================================
base_dir <- "/seurat.obj.comb/cell.anno"
prox_dir <- file.path(base_dir, "/proximity_results")
comb_dir <- file.path(prox_dir, "comb.res")

dir.create(comb_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. Read proximity score files
# =========================================================
read_pscore <- function(subdir, pair_name) {
  df <- read.csv(file.path(prox_dir, subdir, "spatial_pscore.csv"))
  df$cell.pair <- pair_name
  colnames(df)[2] <- "cell.pairs"
  df
}

CAF_FAP_CD8_Tcell.score <- read_pscore("CAF_FAP_CD8_Tcell", "CAF_FAP_CD8_Tcell")
CAF_FAP_CD8_T_GZMB.score <- read_pscore("CAF_FAP_CD8_T_GZMB", "CAF_FAP_CD8_T_GZMB")
Neutrophil_HIF1A_CD8_T_GZMB.score <- read_pscore("Neutrophil_HIF1A_CD8_T_GZMB", "Neutrophil_HIF1A_CD8_T_GZMB")
Neutrophil_HIF1A_CD8_Tcell.score <- read_pscore("Neutrophil_HIF1A_CD8_Tcell", "Neutrophil_HIF1A_CD8_Tcell")
Neutrophil_KI67_CD8_T_GZMB.score <- read_pscore("Neutrophil_KI67_CD8_T_GZMB", "Neutrophil_KI67_CD8_T_GZMB")
Neutrophil_KI67_CD8_Tcell.score <- read_pscore("Neutrophil_KI67_CD8_Tcell", "Neutrophil_KI67_CD8_Tcell")

all.cell.pair.prox.score.cb <- bind_rows(
  CAF_FAP_CD8_Tcell.score,
  CAF_FAP_CD8_T_GZMB.score,
  Neutrophil_HIF1A_CD8_T_GZMB.score,
  Neutrophil_HIF1A_CD8_Tcell.score,
  Neutrophil_KI67_CD8_T_GZMB.score,
  Neutrophil_KI67_CD8_Tcell.score
)

write.csv(
  all.cell.pair.prox.score.cb,
  file = file.path(comb_dir, "all.cell.pair.prox.score.cb.csv"),
  row.names = FALSE
)

# =========================================================
# 3. Prepare annotation columns for plotting
# =========================================================
g.data <- all.cell.pair.prox.score.cb
g.data$spatial_lda_kmeans <- gsub("CN", "H", g.data$spatial_lda_kmeans, fixed = TRUE)
g.data$spatial_lda_kmeans <- factor(g.data$spatial_lda_kmeans, levels = rev(unique(g.data$spatial_lda_kmeans)))

split_cd8 <- function(cell_pair) {
  if (grepl("CD8", cell_pair)) {
    split_result <- sub("_CD8", "|CD8", cell_pair)
    split_result <- unlist(strsplit(split_result, "\\|"))
  } else {
    split_result <- c(cell_pair, NA)
  }
  split_result
}

split_pairs <- t(sapply(g.data$cell.pair, split_cd8))
split_df <- as.data.frame(split_pairs)
colnames(split_df) <- c("CellType1", "CellType2")

g.data <- cbind(g.data, split_df)

# =========================================================
# 4. Plot proximity density heatmap
# =========================================================
min_density <- min(g.data$Proximity.Density, na.rm = TRUE)
max_density <- max(g.data$Proximity.Density, na.rm = TRUE)

p_density <- ggplot(
  g.data,
  aes(
    x = cell.pair,
    y = spatial_lda_kmeans,
    fill = Proximity.Density,
    size = Proximity.Density,
    shape = CellType2
  )
) +
  geom_point(color = "black") +
  scale_shape_manual(values = c("CD8_Tcell" = 22, "CD8_T_GZMB" = 24)) +
  scale_fill_gradientn(
    colors = rev(brewer.pal(5, "RdYlBu")),
    limits = c(min_density, max_density),
    breaks = seq(from = min_density, to = max_density, length.out = 5),
    labels = round(seq(from = min_density, to = max_density, length.out = 5), 2)
  ) +
  scale_size(range = c(1, 10)) +
  labs(
    title = "Proximity Density",
    x = "Cell pair",
    y = "Habitat"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(
  filename = file.path(comb_dir, "Proximity_Density_heatmap.pdf"),
  plot = p_density,
  width = 4.58,
  height = 7.18
)

# =========================================================
# 5. Plot proximity volume heatmap
# =========================================================
min_volume <- min(g.data$Proximity.Volume, na.rm = TRUE)
max_volume <- max(g.data$Proximity.Volume, na.rm = TRUE)

p_volume <- ggplot(
  g.data,
  aes(
    x = cell.pair,
    y = spatial_lda_kmeans,
    fill = Proximity.Volume,
    size = Proximity.Volume,
    shape = CellType2
  )
) +
  geom_point(color = "black") +
  scale_shape_manual(values = c("CD8_Tcell" = 22, "CD8_T_GZMB" = 24)) +
  scale_fill_gradientn(
    colors = rev(brewer.pal(5, "RdYlBu")),
    limits = c(min_volume, max_volume),
    breaks = seq(from = min_volume, to = max_volume, length.out = 5),
    labels = round(seq(from = min_volume, to = max_volume, length.out = 5), 2)
  ) +
  scale_size(range = c(1, 10)) +
  labs(
    title = "Proximity Volume",
    x = "Cell pair",
    y = "Habitat"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(
  filename = file.path(comb_dir, "Proximity_Volume_heatmap.pdf"),
  plot = p_volume,
  width = 4.24,
  height = 7
)