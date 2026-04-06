rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(future)
  library(ggplot2)
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(mclust)
  library(ggrepel)
  library(scRNAtoolVis)
})

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 180 * 1024^3)
options(future.seed = TRUE)
plan("multisession", workers = 10)

# =========================================================
# 1. Define paths and sample-specific QC thresholds
# =========================================================
base_dir <- "/mnt/radonc-li01/private/lyc/CODEX/LUAD"
ws_dir   <- file.path(base_dir, "Whole_slide")
out_dir  <- file.path(base_dir, "results/cell.anno/WS")
qc_dir   <- file.path(out_dir, "qc.post")
comb_dir <- file.path(out_dir, "seurat.obj.comb")
anno_dir <- file.path(out_dir, "cell.annotation")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comb_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(anno_dir, recursive = TRUE, showWarnings = FALSE)


# =========================================================
# 2. Process one whole-slide sample
# =========================================================
process_ws_sample <- function(sample_name, remove_qc = FALSE, marker = NULL, threshold = NULL,
                              ws_dir, qc_dir) {
  message("Processing: ", sample_name)

  in_file <- file.path(ws_dir, sample_name, "quantification", paste0(sample_name, "--unmicst_cellRing.csv"))
  sample_out_dir <- file.path(qc_dir, sample_name)
  dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- fread(in_file, sep = ",", header = TRUE, nThread = 20)

  if (remove_qc && !is.null(marker) && !is.na(marker) && !is.null(threshold) && !is.na(threshold)) {
    del.cell.data <- dat[get(marker) > threshold]

    if (nrow(del.cell.data) > 0) {
      p <- ggplot(del.cell.data, aes(x = X_centroid, y = Y_centroid)) +
        geom_point() +
        scale_y_reverse() +
        theme_classic()
      ggsave(
        filename = file.path(sample_out_dir, paste0(sample_name, "_removed_cells.pdf")),
        plot = p, width = 5, height = 5
      )
      dat <- dat[!CellID %in% del.cell.data$CellID]
    }
  }

  expr_mat <- dat[, c(1, 3:43), with = FALSE]
  meta_mat <- dat[, c(1, 45:53), with = FALSE]

  fwrite(expr_mat, file = file.path(sample_out_dir, "data.exp.matrix.csv"),
         row.names = FALSE, quote = FALSE, sep = ",", nThread = 20)
  fwrite(meta_mat, file = file.path(sample_out_dir, "meta.data.matrix.csv"),
         row.names = FALSE, quote = FALSE, sep = ",", nThread = 20)

  expr_export <- copy(expr_mat)
  meta_export <- copy(meta_mat)
  expr_export$CellID <- paste0(sample_name, "_cell", expr_export$CellID)
  meta_export$CellID <- paste0(sample_name, "_cell", meta_export$CellID)
  meta_export$sample <- sample_name

  fwrite(expr_export, file = file.path(sample_out_dir, "data.exp.matrix.new.csv"),
         row.names = FALSE, quote = FALSE, sep = ",", nThread = 20)
  fwrite(meta_export, file = file.path(sample_out_dir, "meta.data.matrix.new.csv"),
         row.names = FALSE, quote = FALSE, sep = ",", nThread = 20)

  expr_df <- as.data.frame(expr_mat)
  meta_df <- as.data.frame(meta_mat)

  rownames(expr_df) <- paste0(sample_name, "_cell", meta_df$CellID)
  rownames(meta_df) <- paste0(sample_name, "_cell", meta_df$CellID)

  expr_df <- expr_df[, -1, drop = FALSE]
  meta_df <- meta_df[, -1, drop = FALSE]

  colnames(expr_df)[colnames(expr_df) == "Collagen_IV"] <- "Collagen.IV"
  colnames(expr_df)[colnames(expr_df) == "Granzyme_B"]  <- "Granzyme.B"
  colnames(expr_df)[colnames(expr_df) == "E-cadherin"]  <- "E.cadherin"

  count_mat <- t(as.matrix(expr_df))
  seu <- CreateSeuratObject(counts = count_mat, project = "CODEX", meta.data = meta_df)
  seu$sample <- sample_name

  saveRDS(seu, file = file.path(sample_out_dir, "codex.obj.qc.post.rds"))
  return(seu)
}

# =========================================================
# 3. Run sample-level QC and create Seurat objects
# =========================================================
seurat.list <- vector("list", nrow(sample_qc))
names(seurat.list) <- sample_qc$sample

for (i in seq_len(nrow(sample_qc))) {
  seurat.list[[i]] <- process_ws_sample(
    sample_name = sample_qc$sample[i],
    remove_qc   = sample_qc$remove_qc[i],
    marker      = sample_qc$marker[i],
    threshold   = sample_qc$threshold[i],
    ws_dir      = ws_dir,
    qc_dir      = qc_dir
  )
}

# =========================================================
# 4. Merge all whole-slide samples
# =========================================================
combined_seurat <- Reduce(function(x, y) merge(x, y, merge.data = TRUE), seurat.list)
combined_seurat <- JoinLayers(combined_seurat)

if ("sample" %in% colnames(combined_seurat@meta.data)) {
  colnames(combined_seurat@meta.data)[colnames(combined_seurat@meta.data) == "sample"] <- "Sample"
}

saveRDS(combined_seurat, file = file.path(comb_dir, "combined_seurat_object.rds"))

# =========================================================
# 5. Dimensional reduction and clustering
# =========================================================
codex.obj <- combined_seurat

codex.obj <- NormalizeData(
  codex.obj,
  normalization.method = "CLR",
  margin = 2,
  block.size = 1000,
  verbose = TRUE
)

codex.obj <- ScaleData(codex.obj, block.size = 1000, verbose = TRUE)
VariableFeatures(codex.obj) <- rownames(codex.obj)

codex.obj <- RunPCA(
  codex.obj,
  npcs = 20,
  verbose = FALSE,
  seed.use = 1,
  approx = TRUE
)

if (!"batch" %in% colnames(codex.obj@meta.data)) {
  codex.obj$batch <- if ("Sample" %in% colnames(codex.obj@meta.data)) codex.obj$Sample else codex.obj$orig.ident
}

codex.obj <- RunHarmony(
  codex.obj,
  group.by.vars = "batch",
  plot_convergence = FALSE,
  max_iter = 10,
  early_stop = TRUE,
  theta = 4,
  sigma = 0.1,
  ncores = 1
)

codex.obj <- RunUMAP(
  codex.obj,
  dims = 1:20,
  reduction = "harmony",
  umap.method = "uwot",
  n.neighbors = 15,
  min.dist = 0.3,
  metric = "cosine",
  seed.use = 42,
  verbose = TRUE
)

codex.obj <- FindNeighbors(
  codex.obj,
  dims = 1:20,
  reduction = "harmony",
  k.param = 15,
  verbose = TRUE
)

codex.obj <- FindClusters(
  codex.obj,
  resolution = 0.2,
  random.seed = 1,
  algorithm = 2,
  verbose = TRUE
)

pdf(file.path(comb_dir, "dimplot_by_cluster.pdf"), width = 8, height = 6)
print(DimPlot(codex.obj, raster = TRUE, group.by = "seurat_clusters", label = TRUE))
dev.off()

saveRDS(codex.obj, file = file.path(comb_dir, "all.combined.seurat.obj.cluster.post.rds"))

# =========================================================
# 6. Marker analysis for cluster interpretation
# =========================================================
Idents(codex.obj) <- "seurat_clusters"

diff.marker.list <- FindAllMarkers(
  codex.obj,
  only.pos = FALSE,
  logfc.threshold = 0.2
)
diff.marker.list <- diff.marker.list[order(diff.marker.list$avg_log2FC, decreasing = TRUE), ]
write.csv(diff.marker.list, file = file.path(comb_dir, "diff.marker.list.csv"), row.names = FALSE)

pdf(file.path(comb_dir, "featureplot_selected_markers.pdf"), width = 10, height = 8)
print(FeaturePlot(codex.obj, features = c("CD20", "CD8", "PanCK", "CD163", "CD138", "CD11c", "CD68", "CD14", "CD21")))
dev.off()

# =========================================================
# 7. First-pass major cell type annotation
# =========================================================
cell.anno <- as.character(codex.obj$seurat_clusters)
cell.anno[cell.anno == "0"]  <- "Tcell"
cell.anno[cell.anno == "1"]  <- "Malig"
cell.anno[cell.anno == "2"]  <- "Plasma"
cell.anno[cell.anno == "3"]  <- "Macro"
cell.anno[cell.anno == "4"]  <- "Endo"
cell.anno[cell.anno == "5"]  <- "Fibro"
cell.anno[cell.anno == "6"]  <- "Neutro"
cell.anno[cell.anno == "7"]  <- "Malig"
cell.anno[cell.anno == "8"]  <- "Smooth.muscle"
cell.anno[cell.anno == "9"]  <- "Bcell"
cell.anno[cell.anno == "10"] <- "Fibro"
cell.anno[cell.anno == "11"] <- "Others"
cell.anno[cell.anno == "12"] <- "Malig"
cell.anno[cell.anno == "13"] <- "Malig"
cell.anno[cell.anno == "14"] <- "Others"

codex.obj$cell.anno <- cell.anno

pdf(file.path(anno_dir, "dimplot_major_annotation.pdf"), width = 8, height = 6)
print(DimPlot(codex.obj, label = TRUE, raster = TRUE, group.by = "cell.anno", pt.size = 1) + NoLegend())
dev.off()

# =========================================================
# 8. Refine monocyte/macrophage compartment
# =========================================================
Mono.Macro.obj <- subset(codex.obj, subset = seurat_clusters %in% c(3))
VariableFeatures(Mono.Macro.obj) <- rownames(Mono.Macro.obj)

Mono.Macro.obj <- RunPCA(
  Mono.Macro.obj,
  npcs = 20,
  verbose = FALSE,
  seed.use = 1,
  approx = FALSE
)

Mono.Macro.obj <- RunHarmony(
  Mono.Macro.obj,
  group.by.vars = "batch",
  plot_convergence = FALSE,
  ncores = 10
)

Mono.Macro.obj <- FindNeighbors(
  Mono.Macro.obj,
  dims = 1:6,
  reduction = "harmony",
  verbose = FALSE
)

Mono.Macro.obj <- FindClusters(
  Mono.Macro.obj,
  resolution = 0.1,
  algorithm = 1,
  verbose = FALSE
)

Idents(Mono.Macro.obj) <- "seurat_clusters"

mono.marker.list <- FindAllMarkers(
  Mono.Macro.obj,
  only.pos = FALSE,
  logfc.threshold = 0.1
)
mono.marker.list <- mono.marker.list[order(mono.marker.list$avg_log2FC, decreasing = TRUE), ]
write.csv(mono.marker.list, file = file.path(anno_dir, "Mono.Macro.diff.marker.list.csv"), row.names = FALSE)

pdf(file.path(anno_dir, "Mono_Macro_dotplot.pdf"), width = 6, height = 4)
print(DotPlot(Mono.Macro.obj, features = c("CD14", "CD68", "CD163")))
dev.off()

mono.anno <- as.character(Mono.Macro.obj$seurat_clusters)
mono.anno[mono.anno == "0"] <- "Monocyte"
mono.anno[mono.anno == "1"] <- "Macrophage"
Mono.Macro.obj$cell.anno <- mono.anno

codex.obj$cell.anno.refined <- codex.obj$cell.anno
codex.obj$cell.anno.refined[colnames(Mono.Macro.obj)] <- Mono.Macro.obj$cell.anno

saveRDS(Mono.Macro.obj, file = file.path(anno_dir, "Mono.Macro.obj.rds"))

# =========================================================
# 9. Refine plasma-like compartment
# =========================================================
B.plasma.obj <- subset(codex.obj, subset = seurat_clusters %in% c(2))
VariableFeatures(B.plasma.obj) <- rownames(B.plasma.obj)

B.plasma.obj <- RunPCA(
  B.plasma.obj,
  npcs = 20,
  verbose = FALSE,
  seed.use = 1,
  approx = FALSE
)

B.plasma.obj <- RunHarmony(
  B.plasma.obj,
  group.by.vars = "batch",
  plot_convergence = FALSE,
  ncores = 10
)

B.plasma.obj <- FindNeighbors(
  B.plasma.obj,
  dims = 1:3,
  reduction = "harmony",
  verbose = FALSE
)

B.plasma.obj <- FindClusters(
  B.plasma.obj,
  resolution = 0.1,
  algorithm = 1,
  verbose = FALSE
)

Idents(B.plasma.obj) <- "seurat_clusters"

bplasma.marker.list <- FindAllMarkers(
  B.plasma.obj,
  only.pos = FALSE,
  logfc.threshold = 0.1
)
bplasma.marker.list <- bplasma.marker.list[order(bplasma.marker.list$avg_log2FC, decreasing = TRUE), ]
write.csv(bplasma.marker.list, file = file.path(anno_dir, "B.plasma.diff.marker.list.csv"), row.names = FALSE)

pdf(file.path(anno_dir, "B_plasma_dotplot.pdf"), width = 6, height = 4)
print(DotPlot(B.plasma.obj, features = c("CD138", "PanCK", "EpCAM")))
dev.off()

bplasma.anno <- as.character(B.plasma.obj$seurat_clusters)
bplasma.anno[bplasma.anno == "0"] <- "Plasma"
bplasma.anno[bplasma.anno == "1"] <- "Others"
B.plasma.obj$cell.anno <- bplasma.anno

codex.obj$cell.anno.refined[colnames(B.plasma.obj)] <- B.plasma.obj$cell.anno

saveRDS(B.plasma.obj, file = file.path(anno_dir, "B.plasma.obj.rds"))

# =========================================================
# 10. Optional GMM support for monocyte/macrophage separation
# =========================================================
mono.tmp <- subset(Mono.Macro.obj, subset = cell.anno %in% c("Monocyte", "Macrophage"))
cd14_expression <- FetchData(mono.tmp, vars = "CD14")

if (nrow(cd14_expression) > 10) {
  gmm_cd14 <- Mclust(cd14_expression, G = 2)
  mono.tmp$CD14_cluster <- as.factor(gmm_cd14$classification)
  mono.tmp$cell_type_gmm <- ifelse(
    mono.tmp$CD14_cluster == which.max(gmm_cd14$parameters$mean),
    "Monocyte",
    "Macrophage"
  )

  pdf(file.path(anno_dir, "Mono_Macro_GMM_dotplot.pdf"), width = 6, height = 4)
  print(DotPlot(mono.tmp, features = c("CD14", "CD68", "CD163"), group.by = "cell_type_gmm"))
  dev.off()
}

# =========================================================
# 11. Final visualization and export
# =========================================================
codex.obj$cell.anno.final <- codex.obj$cell.anno.refined
Idents(codex.obj) <- "cell.anno.final"

color.value.anno <- c(
  "Malig" = "gray70",
  "Tcell" = "#6BD089",
  "Bcell" = "#3283FE",
  "Plasma" = "#94D8F6",
  "Neutro" = "#CD3278",
  "Endo" = "#C19859",
  "Fibro" = "#FFC000",
  "Monocyte" = "#7A297B",
  "Macrophage" = "#DC45DF",
  "Smooth.muscle" = "#F79646",
  "Others" = "gray90"
)

cluster.order <- c("Malig", "Tcell", "Bcell", "Plasma", "Neutro", "Endo", "Fibro", "Monocyte", "Macrophage", "Smooth.muscle", "Others")
gene.order <- c("PanCK", "E.cadherin", "EpCAM", "CD3e", "Granzyme.B", "CD20", "CD138", "CD66b", "CD31", "FAP", "CD14", "CD163", "aSMA")
color.value.anno.1 <- color.value.anno[cluster.order]

pdf(file.path(anno_dir, "cluster.umap.pdf"), height = 8, width = 7)
print(
  DimPlot(
    codex.obj,
    label = TRUE,
    label.box = TRUE,
    alpha = 0.9,
    raster = TRUE,
    repel = TRUE,
    label.color = "black",
    group.by = "cell.anno.final",
    cols = color.value.anno.1
  ) + NoLegend()
)
dev.off()

pdf(file.path(anno_dir, "mean.heatmap.pdf"), height = 8, width = 7)
print(
  scRNAtoolVis::averageHeatmap(
    object = codex.obj,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    gene.order = gene.order,
    cluster.order = cluster.order,
    annoCol = TRUE,
    myanCol = color.value.anno.1,
    slot = "data",
    group.by = "cell.anno.final"
  )
)
dev.off()

write.csv(
  codex.obj@meta.data,
  file = file.path(anno_dir, "WS_cell_annotation_metadata.csv"),
  row.names = TRUE
)

saveRDS(codex.obj, file = file.path(anno_dir, "WS_cell_annotation_final.rds"))