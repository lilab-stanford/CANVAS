rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(cowplot)
})

# =========================================================
# 1. Define paths
# =========================================================
base_dir <- "/seurat.obj.comb/cell.anno"
cn_dir   <- file.path(base_dir, "TMA.WS.CN.comb.res/results/scimap.CN")
anno_dir <- file.path(base_dir, "each.compartment.anno")
dens_dir <- file.path(base_dir, "cell.density.res")
out_dir  <- file.path(base_dir, "/results")
raw_dir  <- file.path(out_dir, "circlize.rawdata")
km_dir   <- file.path(raw_dir, "each.CN.cell.surv.plot")
sel_dir  <- file.path(raw_dir, "each.CN.cell.surv.plot.select")
comb_dir <- file.path(raw_dir, "select_CN.cell_surv.plot_comb")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(km_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comb_dir, recursive = TRUE, showWarnings = FALSE)

CN.color.value.anno <- c(
  "CN01" = "gray90", "CN02" = "#EFBC20", "CN03" = "#E3F177", "CN04" = "#7CAB7C",
  "CN05" = "#D98D8B", "CN06" = "#F79646", "CN07" = "#8BCCDD", "CN08" = "#4F81BD",
  "CN09" = "#7A297B", "CN10" = "#C785C0"
)

color.value.anno <- c(
  "Cancer.cell" = "#B17A86", "Tcyto" = "#6BD089", "Th" = "#A5CA6B",
  "Treg" = "#D8E063", "NK.like" = "#138535", "Neutrophil" = "#00D5F2",
  "Plasma.cell" = "#75C5D9", "Bcell" = "#2DA2BF", "Mono" = "#7A297B",
  "M1" = "#C785C8", "M2" = "#BDB6D6", "DCs" = "#586B8C",
  "CAF" = "#FFCB00", "Smooth.muscle" = "#E19F73", "Endo" = "#F08B1E",
  "Others" = "gray90"
)

# =========================================================
# 2. Load CN annotation and refined compartment objects
# =========================================================
TMA.WS.CN.results <- fread(
  file.path(cn_dir, "tma.ws.meta.data.CN.anno.kmeans.10.31.final.csv"),
  sep = ","
)

TMA.CN.results <- TMA.WS.CN.results %>%
  filter(str_detect(imageid, fixed("TMA")))
rownames(TMA.CN.results) <- TMA.CN.results$CellID

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

TMA.seurat.obj <- readRDS(file.path(dens_dir, "codex.obj.update.337.rds"))

# =========================================================
# 3. Integrate refined detailed annotations into the TMA object
# =========================================================
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
# 4. Prepare CN metadata and attach CN labels
# =========================================================
CN.metadata <- TMA.CN.results
CN.metadata$spatial_lda_kmeans <- CN.metadata$spatial_lda_kmeans.10

CN.metadata <- as.data.frame(CN.metadata)
colnames(CN.metadata)[colnames(CN.metadata) == "cell.anno"] <- "phenotype"
colnames(CN.metadata)[colnames(CN.metadata) == "imageid"]   <- "spatial.core.id"

CN.metadata <- CN.metadata %>%
  mutate(
    spatial_lda_kmeans = as.integer(str_remove(spatial_lda_kmeans, "CN")),
    spatial_lda_kmeans = spatial_lda_kmeans + 1,
    spatial_lda_kmeans = sprintf("CN%02d", spatial_lda_kmeans)
  )

CN.metadata$spatial_lda_kmeans <- factor(
  as.character(CN.metadata$spatial_lda_kmeans),
  levels = names(CN.color.value.anno)
)

CN.metadata$phenotype <- factor(
  CN.metadata$phenotype,
  levels = unique(names(color.value.anno))
)

CN.metadata <- CN.metadata[, c(
  "CellID", "X_centroid", "Y_centroid",
  "spatial.core.id", "phenotype", "spatial_lda_kmeans"
)]
rownames(CN.metadata) <- CN.metadata$CellID

TMA.seurat.obj <- AddMetaData(TMA.seurat.obj, metadata = CN.metadata, col.name = "spatial_lda_kmeans")
TMA.seurat.metadata <- TMA.seurat.obj@meta.data
TMA.seurat.metadata$spatial_lda_kmeans <- factor(
  TMA.seurat.metadata$spatial_lda_kmeans,
  levels = names(CN.color.value.anno)
)

cell.anno.2st <- TMA.seurat.metadata$cell.anno.2st
cell.anno.2st[cell.anno.2st == "myCAF"] <- "Smooth.muscle"
TMA.seurat.metadata$cell.anno.2st <- cell.anno.2st

fwrite(TMA.seurat.metadata, file.path(out_dir, "TMA.seurat.metadata.csv"), sep = ",", row.names = TRUE)

# =========================================================
# 5. Build CN-by-detailed-cell-type summary matrix
# =========================================================
cell.anno.detail.CN.table.sum.res <- as.data.frame.array(
  table(TMA.seurat.metadata$cell.anno.2st, TMA.seurat.metadata$spatial_lda_kmeans)
)
cell.anno.detail.CN.table.sum.res$cell.anno.2st <- rownames(cell.anno.detail.CN.table.sum.res)

anno_mapping <- unique(TMA.seurat.metadata[, c("cell.anno.2st", "cell.anno.3st")])
write.csv(anno_mapping, file = file.path(raw_dir, "all.cell_anno_mapping.csv"), row.names = FALSE)

combined_data <- merge(
  anno_mapping,
  cell.anno.detail.CN.table.sum.res,
  by = "cell.anno.2st",
  all = TRUE
)
combined_data <- combined_data[order(combined_data$cell.anno.3st), ]

write.csv(combined_data, file = file.path(raw_dir, "cell.anno.detail.CN.table.sum.res.csv"), row.names = FALSE)

cn_data <- as.matrix(combined_data[, 3:12])
column_sums <- colSums(cn_data)
frequency_matrix_summary <- sweep(cn_data + 1, 2, column_sums, FUN = "/")
rownames(frequency_matrix_summary) <- paste(
  combined_data$cell.anno.2st,
  combined_data$cell.anno.3st,
  sep = " | "
)

write.csv(frequency_matrix_summary, file = file.path(raw_dir, "CN_detail_frequency_summary.csv"))

# =========================================================
# 6. Compute overall cell-type abundance and subtype difference
# =========================================================
sam.inf <- read.csv(
  file.path(
    base_dir,
    "cell.density.res_20240812/path.subtype.diff/results/cell.fraction.diff.res/sample.inf.comb.csv"
  ),
  header = TRUE,
  as.is = TRUE
)

sam.inf <- sam.inf[!is.na(sam.inf$type), ]
rownames(sam.inf) <- sam.inf$spatial.core.id
colnames(sam.inf)[colnames(sam.inf) == "spatial_core_id"] <- "spatial.core.id"

total_frequencies <- TMA.seurat.metadata %>%
  count(cell.anno.2st) %>%
  mutate(frequency = (n + 1) / sum(n + 1))

ade_frequencies <- TMA.seurat.metadata %>%
  filter(subtype.y == "adeno") %>%
  count(cell.anno.2st) %>%
  mutate(frequency = (n + 1) / sum(n + 1))

squa_frequencies <- TMA.seurat.metadata %>%
  filter(subtype.y == "squamous cell") %>%
  count(cell.anno.2st) %>%
  mutate(frequency = (n + 1) / sum(n + 1))

colnames(total_frequencies)[2:3] <- c("total.n", "total.frequency")
colnames(ade_frequencies)[2:3]   <- c("ade.n", "ade.frequency")
colnames(squa_frequencies)[2:3]  <- c("squa.n", "squa.frequency")

freq.fc.sum.res <- total_frequencies %>%
  left_join(ade_frequencies, by = "cell.anno.2st") %>%
  left_join(squa_frequencies, by = "cell.anno.2st") %>%
  mutate(
    fold_change = ade.frequency / squa.frequency,
    log2_fold_change = log2(fold_change)
  )

freq.fc.sum.res[is.na(freq.fc.sum.res)] <- 0

cell_counts <- TMA.seurat.metadata %>%
  count(spatial.core.id, cell.anno.2st, subtype.y)

cell_frequencies <- cell_counts %>%
  pivot_wider(names_from = cell.anno.2st, values_from = n, values_fill = list(n = 0))

cell_frequencies_norm <- cell_frequencies %>%
  mutate(across(-c(spatial.core.id, subtype.y), ~ .x / sum(.x, na.rm = TRUE)))

adeno_data <- filter(cell_frequencies_norm, subtype.y == "adeno")
squamous_data <- filter(cell_frequencies_norm, subtype.y == "squamous cell")

p_values <- purrr::map2_df(
  adeno_data[-c(1, 2)],
  squamous_data[-c(1, 2)],
  ~ wilcox.test(.x, .y, exact = TRUE)$p.value
)

p_values_df <- tibble(
  cell_type = names(p_values),
  p_value = unlist(p_values)
)

freq.fc.sum.res <- freq.fc.sum.res %>%
  left_join(p_values_df, by = c("cell.anno.2st" = "cell_type"))

write.csv(freq.fc.sum.res, file = file.path(raw_dir, "freq.fc.sum.res.csv"), row.names = FALSE)

# =========================================================
# 7. Build sample-level CN-by-cell-type count matrix
# =========================================================
frequency_table <- TMA.seurat.metadata %>%
  group_by(spatial.core.id, CN = spatial_lda_kmeans, cell_type = cell.anno.2st) %>%
  summarise(count = n(), .groups = "drop")

frequency_table.raw <- frequency_table

frequency_table_wide <- frequency_table %>%
  unite(col = "CN_cell_type", CN, cell_type, sep = "_") %>%
  select(spatial.core.id, CN_cell_type, count) %>%
  pivot_wider(
    names_from = CN_cell_type,
    values_from = count,
    values_fill = list(count = 0)
  ) %>%
  as.data.frame()

sam.inf.surv <- data.frame(
  spatial.core.id = sam.inf$spatial.core.id,
  status = sam.inf$Survival.status,
  time = sam.inf$surv.time.new.CN.mod
)
sam.inf.surv$status <- ifelse(sam.inf.surv$status == "dead", 1, 0)
rownames(sam.inf.surv) <- sam.inf.surv$spatial.core.id

CN.celltype.freq.mode.data <- merge(
  sam.inf.surv,
  frequency_table_wide,
  by = "spatial.core.id",
  all.x = TRUE
)

write.csv(CN.celltype.freq.mode.data, file = file.path(out_dir, "CN.celltype.freq.mode.data.csv"), row.names = FALSE)

# =========================================================
# 8. Convert counts to within-sample frequencies and filter sparse features
# =========================================================
frequency_matrix <- CN.celltype.freq.mode.data %>%
  dplyr::select(-status, -time) %>%
  mutate(row_total = rowSums(dplyr::select(., -spatial.core.id)))

frequency_matrix <- frequency_matrix %>%
  mutate(across(-c(spatial.core.id, row_total), ~ .x / row_total))

frequency_matrix$row_total <- NULL
frequency_matrix.uni.raw <- frequency_matrix

frequency_matrix_filter <- apply(
  frequency_matrix[, -1],
  2,
  function(x) sum(x > 0) / length(x)
)

freq.filter.id <- names(frequency_matrix_filter[frequency_matrix_filter > 0.1])
frequency_matrix <- frequency_matrix[, c(1, which(colnames(frequency_matrix) %in% freq.filter.id))]

combined_data_surv <- merge(sam.inf, frequency_matrix, by = "spatial.core.id", all.x = TRUE)

surv.data <- cbind.data.frame(
  status = combined_data_surv$Survival.status,
  time = combined_data_surv$surv.time.new.CN.mod,
  combined_data_surv[, grep("CN", colnames(combined_data_surv), fixed = TRUE)]
)

rownames(surv.data) <- combined_data_surv$spatial.core.id
surv.data$status <- ifelse(surv.data$status == "dead", 1, 0)

# =========================================================
# 9. Optimal cutpoint binarization
# =========================================================
surv.opt.cut.data <- matrix(nrow = nrow(surv.data), ncol = ncol(surv.data))

for (i in 3:ncol(surv.data)) {
  sur.cutoff <- surv_cutpoint(
    surv.data,
    time = "time",
    event = "status",
    variables = colnames(surv.data)[i],
    minprop = 0.1
  )
  sur.cutoff <- as.numeric(sur.cutoff$cutpoint[1])

  exp_vec <- surv.data[, i]
  exp_vec[exp_vec < sur.cutoff] <- 0
  exp_vec[exp_vec != 0] <- 1
  surv.opt.cut.data[, i] <- exp_vec
}

colnames(surv.opt.cut.data) <- colnames(surv.data)
rownames(surv.opt.cut.data) <- rownames(surv.data)
surv.opt.cut.data <- as.data.frame(surv.opt.cut.data)
surv.opt.cut.data[, c(1, 2)] <- surv.data[, c(1, 2)]

write.csv(surv.opt.cut.data, file = file.path(out_dir, "surv.opt.cut.data.csv"))

# =========================================================
# 10. Univariate survival analysis
# =========================================================
data <- surv.opt.cut.data

log.rank.p <- c()
values.p <- c()
values.hr <- c()
values.hr_upper <- c()
values.hr_lower <- c()

for (i in 3:ncol(data)) {
  surv.list <- list(
    time = as.numeric(data[, 2]),
    status = data[, 1],
    exp = as.matrix(data[, i])
  )

  km <- survdiff(Surv(time, status) ~ exp, data = surv.list)
  log.value.p <- 1 - pchisq(km[[5]], 1)
  log.rank.p <- c(log.rank.p, log.value.p)

  cox_result <- coxph(Surv(time, status) ~ exp, data = surv.list)
  cox_terms <- summary(cox_result)

  values.p <- c(values.p, signif(cox_terms$wald["pvalue"], digits = 3))
  values.hr <- c(values.hr, signif(cox_terms$coef[2], digits = 3))
  values.hr_upper <- c(values.hr_upper, signif(cox_terms$conf.int[, "upper .95"], digits = 3))
  values.hr_lower <- c(values.hr_lower, signif(cox_terms$conf.int[, "lower .95"], digits = 3))
}

tma.CN.surv.results <- cbind.data.frame(
  CN.cell.id = colnames(data)[-c(1:2)],
  lr.p = log.rank.p,
  uni.cox.p = values.p,
  hr = values.hr,
  hr.upper = values.hr_upper,
  hr.low = values.hr_lower
)

tma.CN.surv.results <- tma.CN.surv.results %>%
  separate(
    CN.cell.id,
    into = c("CN.id", "Cell.id"),
    sep = "_",
    extra = "merge",
    remove = FALSE
  )

write.csv(tma.CN.surv.results, file = file.path(out_dir, "tma.CN.surv.results.csv"), row.names = FALSE)
write.csv(tma.CN.surv.results, file = file.path(raw_dir, "tma.CN.surv.results.csv"), row.names = FALSE)

# =========================================================
# 11. Select representative CEUs for downstream display
# =========================================================
CN01.s.cell.id <- c("CN01_Epi", "CN01_Epi_CD44.GATA3.HLA.A.PDL1+LAG3-", "CN01_Epi_Ecadherin.CD44.HLA.A.HLA.E.HLA.DR+PD.L1.LAG3-", "CN01_Epi_Ki67+")
CN02.s.cell.id <- c("CN02_M1")
CN03.s.cell.id <- c("CN03_aBcell", "CN03_hBcell", "CN03_pBcell")
CN04.s.cell.id <- c("CN04_aCAF", "CN04_mCAF", "CN04_Hipoixa.CAFs")
CN05.s.cell.id <- c("CN05_Plasma_LongAP", "CN05_Plasma_Mature")
CN06.s.cell.id <- c("CN06_Neutrophil_HIF1A", "CN06_Neutrophil_KI67")
CN07.s.cell.id <- c("CN07_Epi_Ki67+", "Epi_Ecadherin.CD44.HLA.A.HLA.E.HLA.DR+PD.L1.LAG3-", "Epi_CD44.HLA.DR+PDL1.LAG3-", "CN07_aCAF", "CN07_AntPresEndo", "CN07_ImmModEndo", "CN07_CD8_cytotoxic_Tcell", "CN07_CD8_simulated_cytotoxic_Tcell", "CD4_Memory_Tcell", "CD4_activated_Tcell", "CN07_fDCs", "CN07_Neutrophil_HIF1A", "CN07_Neutrophil_KI67")
CN08.s.cell.id <- c("CN08_CD4_Memory_Tcell", "CN08_CD4_activated_Tcell")
CN09.s.cell.id <- c("Neutrophil_HLA.DR")
CN10.s.cell.id <- c("CN10_AuxImmEndo", "CN10_ProlImmEndo")

CNs.selec.cell.cb.res <- rbind.data.frame(
  tma.CN.surv.results[tma.CN.surv.results$CN.cell.id %in% CN01.s.cell.id, ],
  tma.CN.surv.results[tma.CN.surv.results$CN.cell.id %in% CN02.s.cell.id, ],
  tma.CN.surv.results[tma.CN.surv.results$CN.cell.id %in% CN03.s.cell.id, ],
  tma.CN.surv.results[tma.CN.surv.results$CN.cell.id %in% CN04.s.cell.id, ],
  tma.CN.surv.results[tma.CN.surv.results$CN.cell.id %in% CN05.s.cell.id, ],
  tma.CN.surv.results[tma.CN.surv.results$CN.cell.id %in% CN06.s.cell.id, ],
  tma.CN.surv.results[tma.CN.surv.results$CN.cell.id %in% CN07.s.cell.id, ],
  tma.CN.surv.results[tma.CN.surv.results$CN.cell.id %in% CN08.s.cell.id, ],
  tma.CN.surv.results[tma.CN.surv.results$CN.cell.id %in% CN09.s.cell.id, ],
  tma.CN.surv.results[tma.CN.surv.results$CN.cell.id %in% CN10.s.cell.id, ]
)

write.csv(CNs.selec.cell.cb.res, file = file.path(raw_dir, "CNs.select.cell.cb.surv.res.csv"), row.names = FALSE)

# =========================================================
# 12. Generate KM plots for all CN-cell features
# =========================================================
safe_filename <- function(x) {
  x <- gsub("[/\\\\:*?\"<>|]", "_", x)
  x <- gsub("\\s+", "_", x)
  x
}

for (i in 3:ncol(data)) {
  col_name <- colnames(data)[i]

  surv_data <- data.frame(
    time = as.numeric(data[, "time"]),
    status = as.numeric(data[, "status"]),
    exp = as.numeric(data[, col_name])
  )

  fit <- survfit(Surv(time, status) ~ exp, data = surv_data)
  km_test <- survdiff(Surv(time, status) ~ exp, data = surv_data)
  log_rank_p <- 1 - pchisq(km_test$chisq, df = 1)

  cox_result <- coxph(Surv(time, status) ~ exp, data = surv_data)
  cox_summary <- summary(cox_result)
  cox_p <- cox_summary$coefficients[1, "Pr(>|z|)"]
  hr <- exp(cox_summary$coefficients[1, "coef"])
  ci_lower <- cox_summary$conf.int[1, "lower .95"]
  ci_upper <- cox_summary$conf.int[1, "upper .95"]

  pdf_filename <- file.path(km_dir, paste0("survival_plot_", safe_filename(col_name), ".pdf"))
  pdf(pdf_filename, width = 7, height = 8.5)

  km_plot <- ggsurvplot(
    fit,
    data = surv_data,
    palette = c("#4F81BD", "#C0504D"),
    risk.table = TRUE,
    tables.height = 0.2,
    tables.theme = theme_classic(),
    pval = formatC(log_rank_p, format = "e", digits = 2),
    conf.int = FALSE,
    xlab = "Time in months",
    ggtheme = theme_classic(),
    ncensor.plot = FALSE,
    pval.coord = c(100, 0.9),
    title = paste("Survival Analysis for", col_name),
    submain = paste(
      "HR[95% CI]:", round(hr, 3), "[", round(ci_lower, 3), "-", round(ci_upper, 3), "];",
      "\nLog-rank P:", formatC(log_rank_p, format = "e", digits = 2),
      "; Cox P:", formatC(cox_p, format = "e", digits = 2)
    )
  )

  print(km_plot)
  dev.off()
}

# =========================================================
# 13. Copy selected KM plots
# =========================================================
files <- list.files(km_dir, full.names = TRUE)

for (cell_id in CNs.selec.cell.cb.res$CN.cell.id) {
  pattern <- paste0(safe_filename(cell_id), ".pdf")
  files_to_copy <- grep(pattern, files, value = TRUE, fixed = TRUE)

  if (length(files_to_copy) > 0) {
    file.copy(files_to_copy, sel_dir, overwrite = TRUE)
  }
}

# =========================================================
# 14. Combine selected KM plots into one panel
# =========================================================
selected.os.id <- c(
  "CN01_Epi", "CN02_M1", "CN03_aBcell", "CN04_Hipoixa.CAFs", "CN05_Plasma_LongAP",
  "CN06_Neutrophil_HIF1A", "CN07_CD8_simulated_cytotoxic_Tcell", "CN07_Epi_Ki67+",
  "CN08_CD4_activated_Tcell", "CN10_AuxImmEndo"
)

pdf(file.path(comb_dir, "combined_survival_plots_10years.pdf"), width = 17.17, height = 9.22)

list_of_plots <- list()

for (col_name in selected.os.id) {
  surv_data <- data.frame(
    time = as.numeric(data[, "time"]),
    status = as.numeric(data[, "status"]),
    exp = as.numeric(data[, col_name])
  )

  fit <- survfit(Surv(time, status) ~ exp, data = surv_data)
  km_test <- survdiff(Surv(time, status) ~ exp, data = surv_data)
  log_rank_p <- 1 - pchisq(km_test$chisq, df = 1)

  cox_result <- coxph(Surv(time, status) ~ exp, data = surv_data)
  cox_summary <- summary(cox_result)
  cox_p <- cox_summary$coefficients[1, "Pr(>|z|)"]
  hr <- exp(cox_summary$coefficients[1, "coef"])
  ci_lower <- cox_summary$conf.int[1, "lower .95"]
  ci_upper <- cox_summary$conf.int[1, "upper .95"]

  km_plot <- ggsurvplot(
    fit,
    data = surv_data,
    palette = c("#4F81BD", "#C0504D"),
    risk.table = TRUE,
    tables.height = 0.2,
    tables.theme = theme_classic(),
    pval = formatC(log_rank_p, format = "e", digits = 2),
    conf.int = TRUE,
    xlab = "Time in months",
    xlim = c(0, 120),
    ggtheme = theme_classic(),
    ncensor.plot = FALSE,
    pval.coord = c(100, 0.9),
    title = paste("Survival Analysis for", col_name),
    submain = paste(
      "HR[95% CI]:", round(hr, 3), "[", round(ci_lower, 3), "-", round(ci_upper, 3), "];",
      "\nLog-rank P:", formatC(log_rank_p, format = "e", digits = 2),
      "; Cox P:", formatC(cox_p, format = "e", digits = 2)
    )
  )

  list_of_plots[[col_name]] <- km_plot$plot
}

combined_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow = 2, ncol = 5)
print(combined_plot)
dev.off()