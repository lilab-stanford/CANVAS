rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(corrplot)
  library(ggplot2)
})

# =========================================================
# 1. Define paths
# =========================================================
base_dir <- "/results"
pred_file <- file.path(base_dir, "TCGA.LUAD.csv")
estimate_file <- file.path(base_dir, "/TCGA.pan.cancer.estimate.res.csv")
cibersortx_file <- file.path(base_dir, "/TCGA_LUAD_CIBERSORTX_full.csv")
xcell_file <- file.path(base_dir, "xcell/xCell_TCGA_RSEM.txt")

out_dir <- file.path(base_dir, "TCGA_frac_cn_comp_res")
plot_dir <- file.path(out_dir, "plot")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. Load TCGA LUAD CN prediction results
# =========================================================
path.meta.final.cb <- fread(pred_file, sep = ",")
path.meta.final.cb$predicted <- as.character(path.meta.final.cb$predicted)

replacement_map <- c(
  "0"  = "CN0",
  "1"  = "CN01",
  "2"  = "CN02",
  "3"  = "CN03",
  "4"  = "CN04",
  "5"  = "CN05",
  "6"  = "CN06",
  "7"  = "CN07",
  "8"  = "CN08",
  "9"  = "CN09",
  "10" = "CN10"
)

path.meta.final.cb$predicted <- replacement_map[path.meta.final.cb$predicted]

# Remove background / excluded class
df <- path.meta.final.cb[path.meta.final.cb$predicted != "CN0", ]
df$tcga_pat_short.id <- substring(df$tcga_pat_id, 1, 12)

# =========================================================
# 3. Summarize CN fractions per patient
# =========================================================
cn_frac <- df %>%
  group_by(tcga_pat_short.id, predicted) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tcga_pat_short.id) %>%
  mutate(frequency = count / sum(count)) %>%
  select(-count) %>%
  pivot_wider(
    names_from = predicted,
    values_from = frequency,
    values_fill = list(frequency = 0)
  ) %>%
  as.data.frame()

cn_cols <- c("CN01", "CN02", "CN03", "CN04", "CN05", "CN06", "CN07", "CN08", "CN09", "CN10")
missing_cn <- setdiff(cn_cols, colnames(cn_frac))
if (length(missing_cn) > 0) {
  for (cn in missing_cn) cn_frac[[cn]] <- 0
}
cn_frac <- cn_frac[, c("tcga_pat_short.id", cn_cols)]

# =========================================================
# 4. Load and merge ESTIMATE results
# =========================================================
tcga.estimate.res <- read.csv(estimate_file, header = TRUE, as.is = TRUE)
tcga.luad.estimate.res <- tcga.estimate.res[tcga.estimate.res$cancer.type == "LUAD", ]
tcga.luad.estimate.res$X <- substring(tcga.luad.estimate.res$X, 1, 12)
tcga.luad.estimate.res$X <- gsub(".", "-", tcga.luad.estimate.res$X, fixed = TRUE)
colnames(tcga.luad.estimate.res)[1] <- "tcga_pat_id"

uniq.sam.id <- intersect(cn_frac$tcga_pat_short.id, tcga.luad.estimate.res$tcga_pat_id)

tcga.cn.pred_est.frac_comb <- cbind.data.frame(
  cn_frac[match(uniq.sam.id, cn_frac$tcga_pat_short.id), ],
  tcga.luad.estimate.res[match(uniq.sam.id, tcga.luad.estimate.res$tcga_pat_id), ]
)

write.csv(
  tcga.cn.pred_est.frac_comb,
  file = file.path(out_dir, "tcga.cn.pred_est.frac_comb.csv"),
  row.names = FALSE
)

# =========================================================
# 5. Load and merge CIBERSORTx results
# =========================================================
tcga.luad_ciberx.abs.frac.res <- read.csv(
  cibersortx_file,
  header = TRUE,
  as.is = TRUE,
  check.names = FALSE
)

inter.sam.id <- intersect(
  tcga.luad_ciberx.abs.frac.res$sample_name,
  tcga.cn.pred_est.frac_comb$tcga_pat_short.id
)

tcga.cn.pred_est.frac_cibx_comb <- cbind.data.frame(
  tcga.cn.pred_est.frac_comb[
    match(inter.sam.id, tcga.cn.pred_est.frac_comb$tcga_pat_short.id),
  ],
  tcga.luad_ciberx.abs.frac.res[
    match(inter.sam.id, tcga.luad_ciberx.abs.frac.res$sample_name),
  ]
)

drop_cols <- c(
  "tcga_pat_id", "cancer.type", "sample_name", "method",
  "P-value", "Correlation", "RMSE", "Absolute score (sig.score)", "ESTIMATEScore"
)
drop_cols <- intersect(drop_cols, colnames(tcga.cn.pred_est.frac_cibx_comb))
if (length(drop_cols) > 0) {
  tcga.cn.pred_est.frac_cibx_comb <- tcga.cn.pred_est.frac_cibx_comb[, !colnames(tcga.cn.pred_est.frac_cibx_comb) %in% drop_cols, drop = FALSE]
}

# =========================================================
# 6. Load and merge xCell results
# =========================================================
tcga_xcell_frac <- read.table(
  xcell_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

tcga_xcell_frac <- t(tcga_xcell_frac)
rownames(tcga_xcell_frac) <- substring(rownames(tcga_xcell_frac), 1, 12)
rownames(tcga_xcell_frac) <- gsub(".", "-", rownames(tcga_xcell_frac), fixed = TRUE)
colnames(tcga_xcell_frac) <- paste0("xc_", colnames(tcga_xcell_frac))

inter.sam.id <- intersect(
  rownames(tcga_xcell_frac),
  tcga.cn.pred_est.frac_cibx_comb$tcga_pat_short.id
)

tcga.cn.pred_est.frac_cibx_xcell_comb <- cbind.data.frame(
  tcga.cn.pred_est.frac_cibx_comb[
    match(inter.sam.id, tcga.cn.pred_est.frac_cibx_comb$tcga_pat_short.id),
  ],
  tcga_xcell_frac[inter.sam.id, , drop = FALSE]
)

if ("tcga_pat_short.id" %in% colnames(tcga.cn.pred_est.frac_cibx_xcell_comb)) {
  cor.data <- tcga.cn.pred_est.frac_cibx_xcell_comb[, !colnames(tcga.cn.pred_est.frac_cibx_xcell_comb) %in% "tcga_pat_short.id", drop = FALSE]
} else {
  cor.data <- tcga.cn.pred_est.frac_cibx_xcell_comb
}

write.csv(
  tcga.cn.pred_est.frac_cibx_xcell_comb,
  file = file.path(out_dir, "tcga.cn.pred_est.frac_cibx_xcell_comb.csv"),
  row.names = FALSE
)

# =========================================================
# 7. Correlation test helper
# =========================================================
cor.mtest <- function(mat, method = "spearman", conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  lowCI.mat <- matrix(NA, n, n)
  uppCI.mat <- matrix(NA, n, n)

  diag(p.mat) <- 0
  diag(lowCI.mat) <- 1
  diag(uppCI.mat) <- 1

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- suppressWarnings(cor.test(mat[, i], mat[, j], method = method, conf.level = conf.level))
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- if (!is.null(tmp$conf.int)) tmp$conf.int[1] else NA
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- if (!is.null(tmp$conf.int)) tmp$conf.int[2] else NA
    }
  }

  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  colnames(lowCI.mat) <- rownames(lowCI.mat) <- colnames(mat)
  colnames(uppCI.mat) <- rownames(uppCI.mat) <- colnames(mat)

  list(p = p.mat, lowCI = lowCI.mat, uppCI = uppCI.mat)
}

# =========================================================
# 8. Compute correlation matrix
# =========================================================
cor.data <- cor.data[, sapply(cor.data, is.numeric), drop = FALSE]

M <- cor(cor.data, method = "spearman", use = "pairwise.complete.obs")
testRes <- cor.mtest(cor.data, method = "spearman", conf.level = 0.95)

# =========================================================
# 9. Plot and export correlation heatmap
# =========================================================
save.id <- file.path(plot_dir, "tcga.luad_est_cib_xcell_frac_cor.plot.pdf")

pdf(save.id, width = 20, height = 20)
corrplot(
  M,
  p.mat = testRes$p,
  method = "color",
  diag = FALSE,
  type = "upper",
  sig.level = c(0.001, 0.01, 0.05),
  pch.cex = 0.9,
  insig = "label_sig",
  pch.col = "grey20",
  order = "original",
  col = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(200)),
  tl.cex = 0.5,
  tl.col = "black"
)
dev.off()

# =========================================================
# 10. Export long-format matrix for optional ggplot use
# =========================================================
g.data <- reshape2::melt(M)
write.csv(g.data, file = file.path(out_dir, "tcga_luad_cn_immune_correlation_long.csv"), row.names = FALSE)