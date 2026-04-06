rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(future)
  library(harmony)
})

options(Seurat.object.assay.version = "v5")
options(future.seed = TRUE)
options(future.globals.maxSize = 180 * 1024^3)
plan("multisession", workers = 20)
set.seed(1)

# =========================
# Paths
# =========================
sample_dir <- "/mnt/radonc-li01/private/lyc/CODEX/LUAD/NSCLC_TMA_2nd/sample.inf"
quant_dir  <- "/mnt/radonc-li01/private/lyc/CODEX/LUAD/TMA_2nd/NSCLC_TMA_Scan1.ome.tiff/quantification"
out_dir    <- "/mnt/radonc-li01/private/lyc/CODEX/LUAD/NSCLC_TMA_2nd/results/cell.anno/anchor.anno"
anno_dir   <- file.path(out_dir, "each.compartment.anno")
seurat_dir <- "/mnt/radonc-li01/private/lyc/CODEX/LUAD/NSCLC_TMA_2nd/results/cell.anno/seurat.creat"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(anno_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(seurat_dir, recursive = TRUE, showWarnings = FALSE)

ref_rds <- "/mnt/radonc-li01/private/lyc/CODEX/LUAD/results/cell.anno/TMA/seurat.obj.comb/cell.anno/cell.density.res/codex.obj.update.337.rds"

# =========================
# 1. Build sample metadata
# =========================
post.inf <- read.csv(
  file.path(sample_dir, "tma_mcmicro_xy_inf.final.csv"),
  header = TRUE, as.is = TRUE
)
post.inf$Pos.id <- paste0(post.inf$Pos.row, post.inf$Pos.col)
post.inf <- post.inf[post.inf$mcmicro.no != "39", ]

sam.inf <- read.csv(
  file.path(sample_dir, "TMA_2nd_sample.inf_20240827.csv"),
  header = TRUE, as.is = TRUE
)
colnames(sam.inf)[1] <- "Pos.id"

sample_map <- merge(post.inf, sam.inf, by = "Pos.id", all = FALSE)
write.csv(
  sample_map,
  file.path(sample_dir, "Scan3_full_centroids.sample.inf.csv"),
  row.names = FALSE
)

sample_map <- sample_map[
  !sample_map$Pathology.diagnosis.new %in%
    c("Lung tissue", "Adjacent normal lung tissue", "Adjacent normal lung tissue (sparse)"),
]

write.csv(
  sample_map,
  file.path(sample_dir, "TMA3_sample.inf.selected.csv"),
  row.names = FALSE
)

# =========================
# 2. Read quantification and create Seurat object
# =========================
file.id <- list.files(quant_dir, pattern = "cellRing", full.names = TRUE)

count.list <- list()
meta.list  <- list()

for (f in file.id) {
  dat <- read.csv(f, header = TRUE, as.is = TRUE)

  core.num <- strsplit(basename(f), "--")[[1]][1]
  meta.data <- dat[, c(1, 45:53)]
  meta.data <- cbind.data.frame(CoreID = core.num, Sample = "core", meta.data)
  meta.data$orig.ident <- "core"
  meta.data$scan.id <- "scan3"
  meta.data$crop.id <- "full"

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

colnames(count.data.comb)[colnames(count.data.comb) == "Collagen_IV"] <- "Collagen.IV"
colnames(count.data.comb)[colnames(count.data.comb) == "Granzyme_B"]  <- "Granzyme.B"

count.mat <- t(as.matrix(count.data.comb))
codex.obj <- CreateSeuratObject(
  counts = count.mat,
  project = "CODEX",
  meta.data = meta.data.comb
)

saveRDS(codex.obj, file.path(seurat_dir, "codex.obj.no.filter.rds"))

# =========================
# 3. Filter selected cores
# =========================
sample.inf.selected <- read.csv(
  file.path(sample_dir, "TMA3_sample.inf.selected.csv"),
  header = TRUE, as.is = TRUE
)

inter.core.id <- intersect(
  unique(sample.inf.selected$mcmicro.no),
  unique(codex.obj@meta.data$CoreID)
)

codex.obj <- subset(codex.obj, subset = CoreID %in% inter.core.id)
saveRDS(codex.obj, file.path(seurat_dir, "codex.obj.scan3.filter.rds"))

# =========================
# 4. Anchor-based annotation
# =========================
ref.codex.obj <- readRDS(ref_rds)

codex.obj <- NormalizeData(codex.obj, normalization.method = "CLR", margin = 2)
codex.obj <- ScaleData(codex.obj)
VariableFeatures(codex.obj) <- rownames(codex.obj)
codex.obj <- RunPCA(codex.obj, npcs = 20, verbose = FALSE, seed.use = 1, approx = FALSE)

# =========================
# 5. Harmony + clustering
# =========================
codex.obj$each.core.id <- paste(codex.obj$scan.id, codex.obj$crop.id, codex.obj$CoreID, sep = "_")

codex.obj <- RunHarmony(
  codex.obj,
  group.by.vars = "each.core.id",
  plot_convergence = FALSE,
  max_iter = 20,
  early_stop = FALSE,
  theta = 4,
  sigma = 0.1,
  ncores = 10
)

codex.obj <- RunUMAP(
  codex.obj,
  dims = 1:20,
  reduction = "harmony",
  verbose = FALSE,
  seed.use = 1
)

codex.obj <- FindNeighbors(
  codex.obj,
  dims = 1:20,
  reduction = "harmony",
  verbose = FALSE
)

codex.obj <- FindClusters(codex.obj, resolution = 0.4, verbose = FALSE, random.seed = 1)
saveRDS(codex.obj, file.path(out_dir, "codex.obj.scan3.cluster.post.rds"))

# =========================
# 6. Major cell annotation
# =========================
cluster_to_major <- c(
  "0"  = "Tumor",
  "1"  = "Mono.Macro",
  "2"  = "Neutrophil",
  "3"  = "Tumor",
  "4"  = "Tcell.CD4",
  "5"  = "Tcell.CD8",
  "6"  = "Endo",
  "7"  = "Fibro",
  "8"  = "Tumor",
  "9"  = "Others",
  "10" = "Others",
  "11" = "Endo",
  "12" = "Tumor",
  "13" = "Plasma",
  "14" = "Tumor",
  "15" = "Tumor",
  "16" = "Bcell",
  "17" = "Neutrophil",
  "18" = "Tumor",
  "19" = "Tumor",
  "20" = "Tumor",
  "21" = "Tumor",
  "22" = "Tumor",
  "23" = "Tumor",
  "24" = "Others",
  "25" = "Others",
  "26" = "Others"
)

codex.obj$cell.anno <- unname(cluster_to_major[as.character(codex.obj$seurat_clusters)])

# =========================
# 7. Smooth muscle refinement
# =========================
Fibro.SM <- subset(codex.obj, subset = seurat_clusters %in% c(7, 13))
aSMA_threshold <- quantile(GetAssayData(Fibro.SM, slot = "data")["aSMA", ], 0.9, na.rm = TRUE)
top10_aSMA_cells <- colnames(Fibro.SM)[GetAssayData(Fibro.SM, slot = "data")["aSMA", ] >= aSMA_threshold]

codex.obj$cell.anno[top10_aSMA_cells] <- "Smooth.muscle"
saveRDS(Fibro.SM, file.path(anno_dir, "Fibro.SM.rds"))

# =========================
# 8. CD4 T-cell refinement
# =========================
Tcell.CD4.obj <- subset(codex.obj, subset = cell.anno == "Tcell.CD4")
Tcell.CD4.obj <- FindNeighbors(Tcell.CD4.obj, dims = 1:20, reduction = "harmony", verbose = FALSE)
Tcell.CD4.obj <- FindClusters(Tcell.CD4.obj, resolution = 0.15, verbose = FALSE, random.seed = 1)

cd4_map <- c("0" = "Th", "1" = "Treg")
Tcell.CD4.obj$cell.anno <- unname(cd4_map[as.character(Tcell.CD4.obj$seurat_clusters)])
codex.obj$cell.anno[colnames(Tcell.CD4.obj)] <- Tcell.CD4.obj$cell.anno

saveRDS(Tcell.CD4.obj, file.path(anno_dir, "Tcell.CD4.obj.rds"))

# =========================
# 9. CD8 / NK refinement
# =========================
Tcell.CD8.obj <- subset(codex.obj, subset = seurat_clusters %in% c(5, 10))
Tcell.CD8.obj <- FindNeighbors(Tcell.CD8.obj, dims = 1:20, reduction = "harmony", verbose = FALSE)
Tcell.CD8.obj <- FindClusters(Tcell.CD8.obj, resolution = 0.3, verbose = FALSE, random.seed = 2)

cd8_map <- c(
  "0" = "Tcyto",
  "1" = "Tcyto",
  "2" = "Others",
  "3" = "Tcyto",
  "4" = "Tcyto",
  "5" = "Nk.like",
  "6" = "Tcyto",
  "7" = "Tcyto",
  "8" = "Tcyto",
  "9" = "Tcyto"
)

Tcell.CD8.obj$cell.anno <- unname(cd8_map[as.character(Tcell.CD8.obj$seurat_clusters)])
codex.obj$cell.anno[colnames(Tcell.CD8.obj)] <- Tcell.CD8.obj$cell.anno

saveRDS(Tcell.CD8.obj, file.path(anno_dir, "Tcell.CD8.obj.rds"))

# =========================
# 10. Monocyte / macrophage refinement
# =========================
mono.macro.obj <- subset(codex.obj, subset = cell.anno == "Mono.Macro")
mono.macro.obj <- FindNeighbors(mono.macro.obj, dims = 1:20, reduction = "harmony", verbose = FALSE)
mono.macro.obj <- FindClusters(mono.macro.obj, resolution = 0.2, verbose = FALSE, random.seed = 1)

mono_map <- c(
  "0" = "M2",
  "1" = "M1",
  "2" = "M2",
  "3" = "Mono",
  "4" = "M1",
  "5" = "DC",
  "6" = "M1",
  "7" = "M1"
)

mono.macro.obj$cell.anno <- unname(mono_map[as.character(mono.macro.obj$seurat_clusters)])
codex.obj$cell.anno[colnames(mono.macro.obj)] <- mono.macro.obj$cell.anno

saveRDS(mono.macro.obj, file.path(anno_dir, "mono.macro.obj.rds"))

# =========================
# 11. Merge clinical metadata
# =========================
TMA3_sample.inf.selected <- read.csv(
  file.path(sample_dir, "TMA3_sample.inf.selected_20240925.csv"),
  header = TRUE, row.names = 1, as.is = TRUE, check.names = FALSE
)

TMA3_sample.inf.selected$each.core.id <- paste0("TMA3_full_", TMA3_sample.inf.selected$mcmicro.no)

TMA3_selected_columns <- TMA3_sample.inf.selected %>%
  select(
    each.core.id, mcmicro.no, Age., Sex., Organ.Anatomic.Site.,
    Pathology.diagnosis., Pathology.diagnosis.new, TNM., Grade.,
    Stage., Type., survival.status., status, time
  )

meta.data <- codex.obj@meta.data
combined_data <- meta.data %>% left_join(TMA3_selected_columns, by = "each.core.id")
rownames(combined_data) <- rownames(meta.data)
codex.obj@meta.data <- combined_data

# =========================
# 12. Standardize labels
# =========================
cell.anno <- as.character(codex.obj$cell.anno)
cell.anno[cell.anno == "Plasma"]  <- "Plasma.cell"
cell.anno[cell.anno == "Nk.like"] <- "NK.like"
cell.anno[cell.anno == "DC"]      <- "DCs"
cell.anno[cell.anno == "Tumor"]   <- "Cancer.cell"
cell.anno[cell.anno == "Fibro"]   <- "CAF"

codex.obj$cell.anno <- cell.anno

# =========================
# 13. Save outputs
# =========================
saveRDS(codex.obj, file.path(anno_dir, "codex.obj.detail.anno.update.rds"))

fwrite(
  codex.obj@meta.data,
  file = file.path(anno_dir, "TMA_2nd_cell.anno.csv"),
  row.names = TRUE, quote = TRUE, sep = ","
)

fwrite(
  codex.obj@meta.data[, c("CoreID", "CellID", "each.core.id", "Pathology.diagnosis.new", "cell.anno.major")],
  file = file.path(anno_dir, "TMA_2nd_major_cell.anno.csv"),
  row.names = TRUE, quote = TRUE, sep = ","
)