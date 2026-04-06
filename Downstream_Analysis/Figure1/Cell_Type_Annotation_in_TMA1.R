rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
})

options(Seurat.object.assay.version = "v5")

# =========================
# Paths
# =========================
base_dir <- "/mnt/radonc-li01/private/lyc/CODEX/LUAD/TMA.20240501"
result_dir <- "/mnt/radonc-li01/private/lyc/CODEX/LUAD/results/cell.anno/TMA.20240613"

dir.create(file.path(result_dir, "scan1.up"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "scan1.low"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "scan2.up"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "scan2.low"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "seurat.obj.comb"), recursive = TRUE, showWarnings = FALSE)

# =========================
# 1. Build sample info tables
# =========================
build_sample_info <- function(sample_file, centroid_file, prefix, out_file) {
  sample.inf <- read.csv(sample_file, header = TRUE, as.is = TRUE)
  centroid.inf <- read.csv(centroid_file, header = TRUE, as.is = TRUE)

  sample.inf$spatial.core.id <- paste0(prefix, "_", sample.inf$x_id, "_", sample.inf$y_id)
  centroid.inf$spatial.core.id <- paste0(prefix, "_", centroid.inf$x_id, "_", centroid.inf$y_id)

  inter.core.id <- intersect(sample.inf$spatial.core.id, centroid.inf$spatial.core.id)

  sample.out <- cbind.data.frame(
    centroid.inf[match(inter.core.id, centroid.inf$spatial.core.id), ],
    sample.inf[match(inter.core.id, sample.inf$spatial.core.id), ]
  )

  write.csv(sample.out, file = out_file, row.names = FALSE)
  sample.out
}

scan1.up.info <- build_sample_info(
  sample_file   = "TMA1_Scan1_crop_up_half.ome.tif/coreograph/sampleinfo_cvt.csv",
  centroid_file = "TMA1_Scan1_crop_up_half.ome.tif/coreograph/centroids_cvt.csv",
  prefix        = "TMA1_up",
  out_file      = file.path(result_dir, "scan1.up", "Scan1_crop_up_centroids.sample.inf.csv")
)

scan1.low.info <- build_sample_info(
  sample_file   = "TMA1_Scan1_crop_low_half.ome.tif/coreograph/sampleinfo_cvt.csv",
  centroid_file = "TMA1_Scan1_crop_low_half.ome.tif/coreograph/centroids_cvt.csv",
  prefix        = "TMA1_low",
  out_file      = file.path(result_dir, "scan1.low", "Scan1_crop_low_centroids.sample.inf.csv")
)

scan2.up.info <- build_sample_info(
  sample_file   = "TMA2_Scan1_crop_up_half.ome.tif/coreograph/sampleinfo_cvt.csv",
  centroid_file = "TMA2_Scan1_crop_up_half.ome.tif/coreograph/centroids_cvt.csv",
  prefix        = "TMA2_up",
  out_file      = file.path(result_dir, "scan2.up", "Scan2_crop_up_centroids.sample.inf.csv")
)

scan2.low.info <- build_sample_info(
  sample_file   = "TMA2_Scan1_crop_low_half.ome.tif/coreograph/sampleinfo_cvt.csv",
  centroid_file = "TMA2_Scan1_crop_low_half.ome.tif/coreograph/centroids_cvt.csv",
  prefix        = "TMA2_low",
  out_file      = file.path(result_dir, "scan2.low", "Scan2_crop_low_centroids.sample.inf.csv")
)

# =========================
# 2. Read one crop and create Seurat object
# =========================
read_crop_as_seurat <- function(quant_dir, scan_id, crop_id, out_dir) {
  file.id <- list.files(quant_dir, pattern = "cellRing", full.names = TRUE)

  count.list <- list()
  meta.list  <- list()

  for (f in file.id) {
    dat <- read.csv(f, header = TRUE, as.is = TRUE)

    core.num <- strsplit(basename(f), "--")[[1]][1]

    meta.data <- dat[, c(1, 45:53)]
    meta.data <- cbind.data.frame(CoreID = core.num, Sample = "core", meta.data)
    meta.data$orig.ident <- "core"
    meta.data$scan.id <- scan_id
    meta.data$crop.id <- crop_id

    rownames(meta.data) <- paste(
      meta.data$scan.id, meta.data$crop.id, meta.data$Sample,
      meta.data$CoreID, paste0("cell", meta.data$CellID), sep = "_"
    )

    count.data <- dat[, 3:43]
    rownames(count.data) <- rownames(meta.data)

    count.list[[core.num]] <- count.data
    meta.list[[core.num]]  <- meta.data
  }

  count.data.comb <- do.call(rbind, count.list)
  meta.data.comb  <- do.call(rbind, meta.list)

  fwrite(count.data.comb,
         file = file.path(out_dir, "count.data.comb.csv"),
         row.names = FALSE, quote = TRUE, sep = ",")

  colnames(count.data.comb)[colnames(count.data.comb) == "Collagen_IV"] <- "Collagen.IV"
  colnames(count.data.comb)[colnames(count.data.comb) == "Granzyme_B"]  <- "Granzyme.B"

  count.mat <- t(as.matrix(count.data.comb))
  seu <- CreateSeuratObject(counts = count.mat, project = "CODEX", meta.data = meta.data.comb)

  saveRDS(seu, file = file.path(out_dir, paste0("codex.obj.", scan_id, ".", crop_id, ".rds")))
  seu
}

codex.obj.scan1.up <- read_crop_as_seurat(
  quant_dir = file.path(base_dir, "TMA1_Scan1_crop_up_half.ome.tiff/quantification"),
  scan_id   = "scan1",
  crop_id   = "up",
  out_dir   = file.path(result_dir, "scan1.up")
)

codex.obj.scan1.low <- read_crop_as_seurat(
  quant_dir = file.path(base_dir, "TMA1_Scan1_crop_low_half.ome.tiff/quantification"),
  scan_id   = "scan1",
  crop_id   = "low",
  out_dir   = file.path(result_dir, "scan1.low")
)

codex.obj.scan2.up <- read_crop_as_seurat(
  quant_dir = file.path(base_dir, "TMA2_Scan1_crop_up_half.ome.tiff/quantification"),
  scan_id   = "scan2",
  crop_id   = "up",
  out_dir   = file.path(result_dir, "scan2.up")
)

codex.obj.scan2.low <- read_crop_as_seurat(
  quant_dir = file.path(base_dir, "TMA2_Scan1_crop_low_half.ome.tiff/quantification"),
  scan_id   = "scan2",
  crop_id   = "low",
  out_dir   = file.path(result_dir, "scan2.low")
)

# =========================
# 3. Add clinical/sample metadata
# =========================
add_sample_metadata <- function(seu, sample_info_file, save_file) {
  clinical.data <- read.csv(sample_info_file, header = TRUE, row.names = 1, as.is = TRUE)
  clinical.data$CoreID <- as.factor(clinical.data$label_id)

  meta.data <- as.data.frame(seu@meta.data)
  merged.meta <- merge(meta.data, clinical.data, by = "CoreID", all.x = TRUE)
  rownames(merged.meta) <- rownames(meta.data)

  seu@meta.data <- merged.meta
  saveRDS(seu, file = save_file)
  seu
}

codex.obj.scan1.up <- add_sample_metadata(
  codex.obj.scan1.up,
  sample_info_file = file.path(result_dir, "scan1.up", "sample.inf.selected.csv"),
  save_file = file.path(result_dir, "scan1.up", "codex.obj.scan1.up.filter.add.meta.rds")
)

codex.obj.scan1.low <- add_sample_metadata(
  codex.obj.scan1.low,
  sample_info_file = file.path(result_dir, "scan1.low", "sample.inf.selected.csv"),
  save_file = file.path(result_dir, "scan1.low", "codex.obj.scan1.low.filter.add.meta.rds")
)

codex.obj.scan2.up <- add_sample_metadata(
  codex.obj.scan2.up,
  sample_info_file = file.path(result_dir, "scan2.up", "sample.inf.selected.csv"),
  save_file = file.path(result_dir, "scan2.up", "codex.obj.scan2.up.filter.add.meta.rds")
)

codex.obj.scan2.low <- add_sample_metadata(
  codex.obj.scan2.low,
  sample_info_file = file.path(result_dir, "scan2.low", "sample.inf.selected.csv"),
  save_file = file.path(result_dir, "scan2.low", "codex.obj.scan2.low.filter.add.meta.rds")
)

# =========================
# 4. Optional subtype cleanup
# =========================
fix_subtype <- function(x) {
  x <- as.character(x)
  x[x == "non-small cell"] <- "adeno"
  x[x == ""] <- "squamous cell"
  x
}

if ("subtype" %in% colnames(codex.obj.scan1.up@meta.data)) {
  codex.obj.scan1.up$subtype  <- fix_subtype(codex.obj.scan1.up$subtype)
  codex.obj.scan1.low$subtype <- fix_subtype(codex.obj.scan1.low$subtype)
  codex.obj.scan2.up$subtype  <- fix_subtype(codex.obj.scan2.up$subtype)
  codex.obj.scan2.low$subtype <- fix_subtype(codex.obj.scan2.low$subtype)
}

saveRDS(codex.obj.scan1.up,  file.path(result_dir, "scan1.up",  "codex.obj.scan1.up.filter.add.meta.rds"))
saveRDS(codex.obj.scan1.low, file.path(result_dir, "scan1.low", "codex.obj.scan1.low.filter.add.meta.rds"))
saveRDS(codex.obj.scan2.up,  file.path(result_dir, "scan2.up",  "codex.obj.scan2.up.filter.add.meta.rds"))
saveRDS(codex.obj.scan2.low, file.path(result_dir, "scan2.low", "codex.obj.scan2.low.filter.add.meta.rds"))

# =========================
# 5. Merge all four crops
# =========================
seurat.objects <- list(
  codex.obj.scan1.up,
  codex.obj.scan1.low,
  codex.obj.scan2.up,
  codex.obj.scan2.low
)

all.combined <- Reduce(function(x, y) merge(x, y), seurat.objects)
all.combined <- JoinLayers(all.combined)

saveRDS(
  all.combined,
  file = file.path(result_dir, "seurat.obj.comb", "all.combined.seurat.obj.rds")
)

meta.data.new <- all.combined@meta.data
write.csv(
  meta.data.new,
  file = file.path(result_dir, "seurat.obj.comb", "meta.data.new.csv"),
  row.names = TRUE
)

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(harmony)
  library(dplyr)
  library(data.table)
  library(mclust)
})

plan("multisession", workers = 20)
options(future.globals.maxSize = 180 * 1024^3)
options(future.seed = TRUE)

setwd("/seurat.obj.comb/cell.anno")

# =========================================================
# 6. Load integrated Seurat object and evaluate clustering
# =========================================================
codex.obj <- readRDS("all.combined.seurat.obj.cluster.post.rds")

# Run an additional clustering resolution for comparison
codex.obj <- FindClusters(
  object = codex.obj,
  verbose = FALSE,
  resolution = 0.5,
  random.seed = 1
)

# Marker detection for two candidate resolutions
Idents(codex.obj) <- "RNA_snn_res.0.2"
markers_res_0.2 <- FindAllMarkers(
  codex.obj,
  only.pos = FALSE,
  logfc.threshold = 0.1
)
markers_res_0.2 <- markers_res_0.2[order(markers_res_0.2$avg_log2FC, decreasing = TRUE), ]
write.csv(markers_res_0.2, file = "RNA_snn_res.0.2.diff.marker.list.csv", row.names = FALSE)

Idents(codex.obj) <- "RNA_snn_res.0.5"
markers_res_0.5 <- FindAllMarkers(
  codex.obj,
  only.pos = FALSE,
  logfc.threshold = 0.1
)
markers_res_0.5 <- markers_res_0.5[order(markers_res_0.5$avg_log2FC, decreasing = TRUE), ]
write.csv(markers_res_0.5, file = "RNA_snn_res.0.5.diff.marker.list.csv", row.names = FALSE)

# Save object with both clustering resolutions
saveRDS(codex.obj, file = "all.combined.seurat.obj.cluster.post.0.2.0.5.rds")

# =========================================================
# 7. First-pass major cell annotation using res.0.2
# =========================================================
codex.obj <- readRDS("all.combined.seurat.obj.cluster.post.0.2.0.5.rds")
Idents(codex.obj) <- "RNA_snn_res.0.2"

cell.anno.1st <- as.character(codex.obj$RNA_snn_res.0.2)

cell.anno.1st[cell.anno.1st == "0"]  <- "Malig"
cell.anno.1st[cell.anno.1st == "1"]  <- "Stromal"
cell.anno.1st[cell.anno.1st == "2"]  <- "Tcell"
cell.anno.1st[cell.anno.1st == "3"]  <- "Endo"
cell.anno.1st[cell.anno.1st == "4"]  <- "T.Bcell"
cell.anno.1st[cell.anno.1st == "5"]  <- "Others"
cell.anno.1st[cell.anno.1st == "6"]  <- "Mono.Macro"
cell.anno.1st[cell.anno.1st == "7"]  <- "Neutro"
cell.anno.1st[cell.anno.1st == "8"]  <- "Tcell"
cell.anno.1st[cell.anno.1st == "9"]  <- "Tcell"
cell.anno.1st[cell.anno.1st == "10"] <- "Plasma"
cell.anno.1st[cell.anno.1st == "11"] <- "Malig"
cell.anno.1st[cell.anno.1st == "12"] <- "Malig"
cell.anno.1st[cell.anno.1st == "13"] <- "Malig"

codex.obj$cell.anno.1st <- cell.anno.1st

# Optional quick marker review for major populations
DotPlot(
  codex.obj,
  features = c(
    "CD14", "CD68", "CD163", "HLA.DR", "CD11c",
    "CD3e", "CD4", "CD8", "FOXP3",
    "CD20", "CD21", "CD138",
    "CD31", "FAP", "aSMA", "MPO", "CD66b",
    "PanCK", "EpCAM", "E.cadherin", "Ki67"
  ),
  group.by = "cell.anno.1st"
)

saveRDS(codex.obj, file = "all.combined.seurat.obj.cell.anno.1st.rds")

# =========================================================
# 8. Refine monocyte/macrophage compartment
# =========================================================
codex.obj <- readRDS("all.combined.seurat.obj.cell.anno.1st.rds")

mono.macro.obj <- subset(codex.obj, subset = cell.anno.1st %in% c("Mono.Macro"))

# Recompute dimensional reduction within the subset
VariableFeatures(mono.macro.obj) <- rownames(mono.macro.obj)

mono.macro.obj <- RunPCA(
  object = mono.macro.obj,
  npcs = 20,
  verbose = FALSE,
  seed.use = 1,
  approx = FALSE
)

mono.macro.obj <- RunHarmony(
  mono.macro.obj,
  group.by.vars = "each.core.id",
  plot_convergence = FALSE,
  iter = 20,
  early_stop = FALSE,
  theta = 4,
  sigma = 0.1,
  ncores = 10
)

mono.macro.obj <- FindNeighbors(
  object = mono.macro.obj,
  dims = 1:20,
  reduction = "harmony",
  verbose = FALSE,
  features = c("CD14", "CD68", "CD163", "HLA.DR", "CD11c")
)

mono.macro.obj <- FindClusters(
  object = mono.macro.obj,
  verbose = FALSE,
  resolution = 0.05,
  random.seed = 1,
  algorithm = 1
)

Idents(mono.macro.obj) <- "RNA_snn_res.0.05"

DotPlot(
  mono.macro.obj,
  features = c("CD14", "CD68", "CD163", "HLA.DR", "CD11c")
)

mono.macro.anno <- as.character(mono.macro.obj$RNA_snn_res.0.05)
mono.macro.anno[mono.macro.anno == "0"] <- "M1"
mono.macro.anno[mono.macro.anno == "1"] <- "Mono"
mono.macro.anno[mono.macro.anno == "2"] <- "M2"

mono.macro.obj$mono.macro.anno <- mono.macro.anno

saveRDS(mono.macro.obj, file = "mono.macro.obj.rds")

# =========================================================
# 9. Optional GMM-based support for Mono vs M2 separation
# =========================================================
mono.tmp <- subset(mono.macro.obj, subset = mono.macro.anno.refined %in% c("Mono", "M2"))

cd14_expression <- FetchData(mono.tmp, vars = "CD14")
gmm_cd14 <- Mclust(cd14_expression, G = 2)

mono.tmp$CD14_cluster <- as.factor(gmm_cd14$classification)
mono.tmp$cell_type_gmm <- ifelse(
  mono.tmp$CD14_cluster == which.max(gmm_cd14$parameters$mean),
  "Mono",
  "M2"
)

DotPlot(
  mono.tmp,
  features = c("CD14", "CD68", "CD163"),
  group.by = "cell_type_gmm"
)

# =========================================================
# 10. Integrate refined mono/macrophage labels back to all cells
# =========================================================
codex.obj$cell.anno.2nd <- codex.obj$cell.anno.1st
codex.obj$cell.anno.2nd[colnames(mono.macro.obj)] <- mono.macro.obj$cell.anno.1st

Idents(codex.obj) <- "cell.anno.2nd"

DotPlot(
  codex.obj,
  features = c(
    "CD14", "CD68", "CD163", "HLA.DR", "CD11c",
    "CD3e", "CD4", "CD8", "FOXP3",
    "CD20", "CD21", "CD138",
    "CD31", "FAP", "aSMA", "MPO", "CD66b",
    "PanCK", "EpCAM", "E.cadherin", "Ki67"
  ),
  group.by = "cell.anno.2nd"
)

saveRDS(codex.obj, file = "all.combined.seurat.obj.cell.anno.refined.rds")

# =========================================================
# 11. Export metadata
# =========================================================
write.csv(
  codex.obj@meta.data,
  file = "all.combined.seurat.obj.cell.anno.refined.metadata.csv",
  row.names = TRUE
)