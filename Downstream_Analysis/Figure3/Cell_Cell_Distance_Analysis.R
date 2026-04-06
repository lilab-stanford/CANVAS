rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

# =========================================================
# 1. Define paths
# =========================================================
base_dir <- "/seurat.obj.comb/cell.anno"
anno_dir <- file.path(base_dir, "each.compartment.anno")
cn_dir   <- file.path(base_dir, "TMA.WS.CN.comb.res/results/results/sumplot")
dens_dir <- file.path(base_dir, "cell.density.res")
prox_dir <- file.path(base_dir, "/CN_details.cell_proximity")
raw_dir  <- file.path(prox_dir, "rawdata")
res_dir  <- file.path(prox_dir, "results")
plot_dir <- file.path(prox_dir, "plot")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. Load detailed annotation objects
# =========================================================
mono.macro.obj <- readRDS(file.path(anno_dir, "mono.macro.obj.rds"))
fibro.obj      <- readRDS(file.path(anno_dir, "fibro.obj.rds"))
endo.obj       <- readRDS(file.path(anno_dir, "endo.obj.rds"))
B.obj          <- readRDS(file.path(anno_dir, "B.obj.rds"))
Plasma.obj     <- readRDS(file.path(anno_dir, "Plasma.obj.rds"))
Neutrophil.obj <- readRDS(file.path(anno_dir, "Neutrophil.obj.rds"))
Tcell.obj      <- readRDS(file.path(anno_dir, "Tcell.obj.rds"))
Prof.cell.obj  <- readRDS(file.path(anno_dir, "Prof.cell.obj.rds"))
Tumor.cell.obj <- readRDS(file.path(anno_dir, "Tumor.cell.obj.rds"))
DC.cell.obj    <- readRDS(file.path(anno_dir, "DC.cell.obj.rds"))

Tcell.obj@meta.data$cell.anno.2st   <- Tcell.obj@meta.data$Tcell.anno
DC.cell.obj@meta.data$cell.anno.2st <- DC.cell.obj@meta.data$DC.cell.anno

# =========================================================
# 3. Load main TMA object and add detailed annotations
# =========================================================
TMA.seurat.obj <- readRDS(file.path(dens_dir, "codex.obj.update.337.rds"))

TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = mono.macro.obj@meta.data, col.name = "cell.anno.2st")
TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = fibro.obj@meta.data,      col.name = "cell.anno.2st")
TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = endo.obj@meta.data,       col.name = "cell.anno.2st")
TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = B.obj@meta.data,          col.name = "cell.anno.2st")
TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = Plasma.obj@meta.data,     col.name = "cell.anno.2st")
TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = Neutrophil.obj@meta.data, col.name = "cell.anno.2st")
TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = Tcell.obj@meta.data,      col.name = "cell.anno.2st")
TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = Prof.cell.obj@meta.data,  col.name = "cell.anno.2st")
TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = Tumor.cell.obj@meta.data, col.name = "cell.anno.2st")
TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = DC.cell.obj@meta.data,    col.name = "cell.anno.2st")

# =========================================================
# 4. Add CN labels
# =========================================================
TMA.CN.results <- fread(file.path(cn_dir, "CN.metadata.csv"), sep = ",", header = TRUE, nThread = 20)
TMA.CN.results <- as.data.frame(TMA.CN.results)
rownames(TMA.CN.results) <- TMA.CN.results$CellID
TMA.CN.results <- TMA.CN.results[, c("CellID", "spatial_lda_kmeans")]

TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = TMA.CN.results, col.name = "spatial_lda_kmeans")

# =========================================================
# 5. Export metadata and counts for downstream distance analysis
# =========================================================
saveRDS(TMA.seurat.obj, file = file.path(raw_dir, "TMA.seurat.obj_detail.cells_CNs.rds"))

comb.metadata <- TMA.seurat.obj@meta.data
fwrite(comb.metadata, file = file.path(raw_dir, "TMA.seurat.obj_detail.cells_CNs_metadata.csv"),
       row.names = TRUE, quote = TRUE, sep = ",")

tma.count.data <- as.data.frame(t(TMA.seurat.obj@assays$RNA$counts))
tma.count.data <- cbind.data.frame(cell.id.uniq = rownames(tma.count.data), tma.count.data)
fwrite(tma.count.data, file = file.path(raw_dir, "tma.count.data.csv"),
       row.names = TRUE, quote = TRUE, sep = ",")

# =========================================================
# 6. Load spatial distance matrix
# =========================================================
spatial_distance <- fread(
  file.path(res_dir, "spatial_distance.csv"),
  sep = ",", header = TRUE, na.strings = "NA"
)
spatial_distance <- as.data.frame(spatial_distance)
rownames(spatial_distance) <- spatial_distance$cell.id.uniq
spatial_distance <- spatial_distance[, -1]
spatial_distance <- log2(spatial_distance + 1)

inter.cell.id <- intersect(rownames(comb.metadata), rownames(spatial_distance))
meta_distance <- cbind.data.frame(
  comb.metadata[inter.cell.id, ],
  spatial_distance[inter.cell.id, ]
)

fwrite(meta_distance, file = file.path(res_dir, "tma1.meta.spatdis.res.csv"),
       row.names = TRUE, quote = TRUE, sep = ",")

detail.cell.id <- colnames(spatial_distance)

# =========================================================
# 7. Density plots for selected cell-pair distance distributions
# =========================================================
my_colors <- c(
  "CN01" = "#E5E5E5",
  "CN07" = "#8BCCDD",
  "Others" = "#ECBC27"
)

from.cell.ids <- c(
  "CAF_HIF1A", "CAF_FAP", "Neutrophil_KI67",
  "Neutrophil_HIF1A", "Tumor_PanCK",
  "Tumor_CD44+PDL1+LAG3-", "Tumor_Ki67"
)
to.cell.ids <- c("CD8_Tcell", "CD8_T_GZMB")
cn.keep <- c("CN01", "CN07")

for (from.id in from.cell.ids) {
  for (to.id in to.cell.ids) {
    g.data <- meta_distance[meta_distance$raw.cell.anno.id == from.id, ]
    g.data$CN.new.label <- ifelse(g.data$spatial_lda_kmeans %in% cn.keep, g.data$spatial_lda_kmeans, "Others")

    if (nrow(g.data) > 0 && to.id %in% colnames(g.data) && !all(is.na(g.data[[to.id]]))) {
      p <- ggplot(g.data, aes(x = !!sym(to.id), fill = CN.new.label, color = CN.new.label)) +
        geom_density(alpha = 0) +
        scale_color_manual(values = my_colors) +
        scale_fill_manual(values = my_colors) +
        labs(title = paste0(from.id, " >> ", to.id), x = to.id, y = "Density") +
        theme_classic()

      ggsave(
        filename = file.path(plot_dir, paste0("DistanceDensity_", from.id, "_to_", to.id, ".pdf")),
        plot = p,
        width = 4.97,
        height = 2.43
      )
    }
  }
}

# =========================================================
# 8. Compare CN01 vs CN07 mean distance for all cell pairs
# =========================================================
data <- meta_distance[, c("raw.cell.anno.id", "cell.anno.3st", "spatial_lda_kmeans", "spatial.core.id", detail.cell.id)]
data <- data[data$spatial_lda_kmeans %in% c("CN01", "CN07"), ]

data.long <- reshape2::melt(data, na.rm = FALSE, value.name = "distance")
colnames(data.long) <- c("from.cell.id", "major.cell.id", "CN.id", "spatial.core.id", "to.cell.id", "distance")
data.long$from.to.cell.id <- paste0(data.long$from.cell.id, ">", data.long$to.cell.id)

average_distances <- data.long %>%
  group_by(CN.id, from.to.cell.id, spatial.core.id) %>%
  summarise(mean_distance = mean(distance, na.rm = TRUE), .groups = "drop")

distances_cn01 <- average_distances %>% filter(CN.id == "CN01")
distances_cn07 <- average_distances %>% filter(CN.id == "CN07")

common_pairs <- intersect(distances_cn01$from.to.cell.id, distances_cn07$from.to.cell.id)

distances_cn01 <- as.data.frame(distances_cn01 %>% filter(from.to.cell.id %in% common_pairs))
distances_cn07 <- as.data.frame(distances_cn07 %>% filter(from.to.cell.id %in% common_pairs))

fc <- numeric(length(common_pairs))
pvalue <- numeric(length(common_pairs))

for (i in seq_along(common_pairs)) {
  cn01_dis <- distances_cn01 %>% filter(from.to.cell.id == common_pairs[i]) %>% pull(mean_distance)
  cn07_dis <- distances_cn07 %>% filter(from.to.cell.id == common_pairs[i]) %>% pull(mean_distance)

  if (length(na.omit(cn01_dis)) < 3 || length(na.omit(cn07_dis)) < 3) {
    fc[i] <- NA
    pvalue[i] <- NA
  } else {
    fc[i] <- mean(cn07_dis, na.rm = TRUE) - mean(cn01_dis, na.rm = TRUE)
    pvalue[i] <- tryCatch(
      wilcox.test(cn01_dis, cn07_dis, exact = FALSE, correct = TRUE)$p.value,
      error = function(e) NA
    )
  }
}

result <- data.frame(
  from_to_cell_id = common_pairs,
  fold_change = fc,
  p_value = pvalue
)
result$fdr <- p.adjust(result$p_value, method = "BH")

write.csv(result, file = file.path(res_dir, "CN01.CN07_detail.cell_diff.csv"), row.names = FALSE)

# =========================================================
# 9. Volcano plot helper for selected source cells
# =========================================================
plot_distance_volcano <- function(result_df, from_pattern, labels_to_keep, out_prefix) {
  plot_df <- result_df[grep(paste0(from_pattern, ">"), result_df$from_to_cell_id, fixed = TRUE), ]
  plot_df <- plot_df[order(plot_df$fold_change, decreasing = TRUE), ]
  plot_df <- plot_df %>%
    separate(from_to_cell_id, into = c("from_cell", "to_cell"), sep = ">", remove = FALSE)

  plot_df <- na.omit(plot_df)
  plot_df$lg.padj <- -log10(plot_df$fdr)
  plot_df$log2.fc <- -plot_df$fold_change

  label_df <- plot_df[plot_df$from_to_cell_id %in% labels_to_keep, ]

  p <- ggplot(plot_df, aes(x = log2.fc, y = lg.padj)) +
    geom_point(aes(color = log2.fc, size = lg.padj), alpha = 0.75) +
    scale_color_gradientn(colors = rev(brewer.pal(5, "RdYlBu"))) +
    scale_size_continuous(range = c(0.1, 3)) +
    geom_text_repel(
      data = label_df,
      aes(label = to_cell),
      color = "red",
      size = 3,
      nudge_x = 0.5,
      segment.color = "grey30",
      segment.size = 0.4
    ) +
    geom_point(data = label_df, aes(x = log2.fc, y = lg.padj, size = lg.padj),
               shape = 21, color = "black") +
    theme_classic() +
    labs(x = "Distance difference", y = "-Log10 adjusted P-value") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme(legend.position = "right", legend.title = element_blank())

  out_subdir <- file.path(res_dir, out_prefix)
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)

  ggsave(
    filename = file.path(out_subdir, paste0(out_prefix, "_distance_volcano.pdf")),
    plot = p,
    width = 5.48,
    height = 4.66
  )
}

plot_distance_volcano(
  result,
  from_pattern = "CAF_FAP",
  labels_to_keep = c(
    "CAF_FAP>Tumor_CD44+PDL1.LAG3-",
    "CAF_FAP>Neutrophil_KI67",
    "CAF_FAP>Tumor_Ki67",
    "CAF_FAP>CD8_T_GZMB",
    "CAF_FAP>CD8_T_KI67",
    "CAF_FAP>CD4_Tcell",
    "CAF_FAP>CD8_Tcell",
    "CAF_FAP>Bcell_HLA.A&E",
    "CAF_FAP>Bcell_HIF1A"
  ),
  out_prefix = "CAF_FAP_CN01_CN07"
)

plot_distance_volcano(
  result,
  from_pattern = "Neutrophil_HIF1A",
  labels_to_keep = c(
    "Neutrophil_HIF1A>Bcell_HLA.A&E",
    "Neutrophil_HIF1A>cDCs",
    "Neutrophil_HIF1A>CD4_Tcell",
    "Neutrophil_HIF1A>CD4_T_KI67",
    "Neutrophil_HIF1A>CD8_Tcell",
    "Neutrophil_HIF1A>CD4_T_CD45RO.PD.1",
    "Neutrophil_HIF1A>CD8_T_GZMB",
    "Neutrophil_HIF1A>Tumor_PanCK",
    "Neutrophil_HIF1A>Tumor_CD44+PDL1+LAG3-",
    "Neutrophil_HIF1A>Tumor_BCL2"
  ),
  out_prefix = "Neutrophil_HIF1A_CN01_CN07"
)

plot_distance_volcano(
  result,
  from_pattern = "Neutrophil_KI67",
  labels_to_keep = c(
    "Neutrophil_KI67>Bcell_HLA.A&E",
    "Neutrophil_KI67>Endo_HLA.DR",
    "Neutrophil_KI67>cDCs",
    "Neutrophil_KI67>CD4_Tcell",
    "Neutrophil_KI67>CD8_Tcell",
    "Neutrophil_KI67>CD4_T_CD45RO.PD.1",
    "Neutrophil_KI67>CD8_T_GZMB",
    "Neutrophil_KI67>Tumor_CD44+PDL1+LAG3-",
    "Neutrophil_KI67>Tumor_HIF1A.FOXP3",
    "Neutrophil_KI67>Tumor_Ki67"
  ),
  out_prefix = "Neutrophil_KI67_CN01_CN07"
)