rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(purrr)
  library(data.table)
  library(ggplot2)
  library(effsize)
  library(cowplot)
})

# =========================================================
# 1. Define paths
# =========================================================
base_dir <- "/cell.anno/cell.density.res"
subtype_dir <- file.path(base_dir, "histological.subtype")

detail_dir <- "/Para_bubble_plot/results"
pair_dir   <- "/Others_pair_comp/results"
final_dir  <- "/Others_pair_comp/results_final"
plot_dir   <- "/Others_pair_comp/plot"

dir.create(pair_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)

eps <- 1e-12

# =========================================================
# 2. Load Seurat object and sample metadata
# =========================================================
TMA.codex.obj <- readRDS(file.path(base_dir, "codex.obj.update.337.rds"))

sample.inf <- read.csv(
  file = file.path(subtype_dir, "sample.inf.comb_20240701.new_his.subtype_surv.mode.csv"),
  header = TRUE,
  row.names = 1,
  as.is = TRUE
)

sample.inf <- sample.inf[!is.na(sample.inf$subtype), ]
sample.inf[sample.inf$subtype == "squamous cell", "histological.subtype"] <-
  sample.inf[sample.inf$subtype == "squamous cell", "histological.subtype.pre"]

drop_cols <- intersect(colnames(sample.inf), c("EGFR", "ALK", "ROS1", "KRAS"))
if (length(drop_cols) > 0) {
  sample.inf <- sample.inf[, !colnames(sample.inf) %in% drop_cols, drop = FALSE]
}

mut.info <- read.csv(
  file = file.path(subtype_dir, "tma_mut_tum.size_uniq.id_comb_20241205.csv"),
  header = TRUE,
  row.names = 1,
  as.is = TRUE
)

mut.info <- mut.info[, c("MRN", "Donor.ID", "tumor.size", "Mut.dect", "Gene.mut")]
sample.inf <- merge(sample.inf, mut.info, by.x = "Donor.no.", by.y = "Donor.ID", all = TRUE)

# =========================================================
# 3. Build sample-level metadata table
# =========================================================
meta.data <- TMA.codex.obj@meta.data
meta.data$study_id <- gsub(" ", "_", meta.data$study_id, fixed = TRUE)

meta.unique <- meta.data %>%
  select(
    spatial.core.id, core_id, study_id, stt_num, tissue, type, subtype,
    spatial.core.id.1, each.core.id
  ) %>%
  distinct(study_id, .keep_all = TRUE)

rownames(meta.unique) <- meta.unique$spatial.core.id
colnames(meta.unique)[colnames(meta.unique) == "study_id"] <- "Study.ID"

sample.unique <- sample.inf %>%
  filter(!is.na(tissue)) %>%
  select(
    spatial.core.id, core_id, Study.ID, stt_num, tissue, type, subtype,
    histological.subtype, spatial.core.id.1, each.core.id,
    Age, Gender, Survival.status, surv.time.new.CN.mod,
    Grade.new, Stage, Invasion, MCT.status, MCT.hilo,
    tumor.size, Mut.dect, Gene.mut
  ) %>%
  distinct(Study.ID, .keep_all = TRUE)

rownames(sample.unique) <- sample.unique$spatial.core.id
sample.inf.comb <- merge(meta.unique, sample.unique, by = "Study.ID", all = TRUE)

# =========================================================
# 4. Compute core-level cell fractions
# =========================================================
counts <- meta.data %>%
  group_by(spatial.core.id.1, cell.anno.3st) %>%
  summarise(count = n(), .groups = "drop")

totals <- counts %>%
  group_by(spatial.core.id.1) %>%
  summarise(total = sum(count), .groups = "drop")

proportions <- counts %>%
  left_join(totals, by = "spatial.core.id.1") %>%
  mutate(proportion = count / total)

wide_data <- proportions %>%
  select(-proportion) %>%
  pivot_wider(
    names_from = cell.anno.3st,
    values_from = count,
    values_fill = list(count = 0)
  ) %>%
  rename(spatial_core_id = spatial.core.id.1)

wide_data <- wide_data %>%
  mutate(across(-c(spatial_core_id, total), ~ .x / total * 100)) %>%
  select(-total)

sample.inf.final <- sample.inf %>%
  rename(spatial_core_id = spatial.core.id) %>%
  filter(!is.na(tissue))

merged_data <- left_join(sample.inf.final, wide_data, by = "spatial_core_id")

detail_cell_types <- c(
  "Cancer.cell", "Bcell", "CAF", "Smooth.muscle", "DCs", "Endo",
  "M1", "M2", "Mono", "NK.like", "Neutrophil", "Others",
  "Plasma.cell", "Tcyto", "Th", "Treg"
)

long_data <- melt(
  merged_data,
  id.vars = c(
    "Study.ID", "spatial_core_id", "core_id", "stt_num", "tissue", "type",
    "subtype", "histological.subtype", "spatial.core.id.1", "Gender",
    "Grade.new", "Survival.status", "surv.time.new.CN.mod", "MCT.count",
    "MCT.status", "MCT.hilo", "neg", "MCT.low..1.9.", "tumor.size",
    "Mut.dect", "Gene.mut"
  ),
  measure.vars = intersect(detail_cell_types, colnames(merged_data)),
  variable.name = "cell_type",
  value.name = "proportion"
)

# =========================================================
# 5. Pairwise comparison across all histological subtypes
# =========================================================
subtype_combinations <- combn(unique(long_data$histological.subtype), 2, simplify = FALSE)

pairwise_subtype_comparison_all_cells <- map_dfr(subtype_combinations, function(pair) {
  s1 <- pair[1]
  s2 <- pair[2]

  long_data %>%
    filter(histological.subtype %in% c(s1, s2)) %>%
    group_by(cell_type) %>%
    reframe({
      x <- proportion[histological.subtype == s1]
      y <- proportion[histological.subtype == s2]

      wt <- tryCatch(
        wilcox.test(x, y, exact = FALSE, correct = FALSE),
        error = function(e) NULL
      )

      W  <- if (is.null(wt)) NA_real_ else unname(wt$statistic)
      pv <- if (is.null(wt)) NA_real_ else wt$p.value

      cd <- suppressWarnings(tryCatch(
        as.numeric(effsize::cliff.delta(x, y, method = "asymptotic")$estimate),
        error = function(e) NA_real_
      ))

      fc <- (mean(x, na.rm = TRUE) + eps) / (mean(y, na.rm = TRUE) + eps)

      tibble(
        Fold_Change = fc,
        PValue = pv,
        W = W,
        CliffDelta = cd,
        group1 = s1,
        group2 = s2
      )
    }) %>%
    mutate(Comparison = paste0(group2, "_vs_", group1)) %>%
    select(cell_type, Fold_Change, PValue, W, CliffDelta, Comparison)
}) %>%
  rename(Cell_Type = cell_type) %>%
  mutate(fdr = p.adjust(PValue, method = "fdr"))

write.csv(
  pairwise_subtype_comparison_all_cells,
  file = file.path(pair_dir, "pairwise_subtype_comparison_all_cells.csv"),
  row.names = FALSE
)

# =========================================================
# 6. LUAD vs LUSC comparison for detailed cell types
# =========================================================
cmp_g1 <- "adeno"
cmp_g2 <- "squamous cell"

luad_lusc_cells <- long_data %>%
  filter(subtype %in% c(cmp_g1, cmp_g2)) %>%
  mutate(subtype = factor(subtype, levels = c(cmp_g1, cmp_g2)))

LUAD_vs_LUSC_celltype_comparison <- luad_lusc_cells %>%
  group_by(cell_type) %>%
  summarise(
    W = tryCatch(
      unname(wilcox.test(proportion ~ subtype, exact = FALSE, correct = FALSE)$statistic),
      error = function(e) NA_real_
    ),
    PValue = tryCatch(
      wilcox.test(proportion ~ subtype, exact = FALSE, correct = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    mean_g1 = mean(proportion[subtype == cmp_g1], na.rm = TRUE),
    mean_g2 = mean(proportion[subtype == cmp_g2], na.rm = TRUE),
    Fold_Change = (mean_g2 + eps) / (mean_g1 + eps),
    CliffDelta = suppressWarnings(tryCatch(
      as.numeric(
        effsize::cliff.delta(
          proportion[subtype == cmp_g2],
          proportion[subtype == cmp_g1],
          method = "asymptotic"
        )$estimate
      ),
      error = function(e) NA_real_
    )),
    .groups = "drop"
  ) %>%
  mutate(
    fdr = p.adjust(PValue, method = "fdr"),
    Comparison = paste0(cmp_g1, "_vs_", cmp_g2)
  ) %>%
  transmute(
    Cell_Type = factor(cell_type),
    Fold_Change,
    PValue,
    W,
    CliffDelta,
    Comparison,
    fdr
  )

write.csv(
  LUAD_vs_LUSC_celltype_comparison,
  file = file.path(pair_dir, "LUAD_vs_LUSC_celltype_comparison.csv"),
  row.names = FALSE
)

# =========================================================
# 7. Define compartment-level summaries
# =========================================================
cell_categories <- list(
  Lymphoid = c("Th", "Tcyto", "Treg", "Bcell", "Plasma.cell", "NK.like"),
  Myeloid  = c("Neutrophil", "DCs", "M2", "M1", "Mono"),
  Stromal  = c("Endo", "CAF", "Smooth.muscle"),
  Immune   = c("Th", "Tcyto", "Treg", "Bcell", "Plasma.cell", "NK.like",
               "Neutrophil", "DCs", "M2", "M1", "Mono")
)

wide_data_comp <- wide_data %>%
  mutate(
    Lymphoid = rowSums(select(., all_of(cell_categories$Lymphoid)), na.rm = TRUE),
    Myeloid  = rowSums(select(., all_of(cell_categories$Myeloid)), na.rm = TRUE),
    Stromal  = rowSums(select(., all_of(cell_categories$Stromal)), na.rm = TRUE),
    Immune   = rowSums(select(., all_of(cell_categories$Immune)), na.rm = TRUE)
  ) %>%
  select(spatial_core_id, Lymphoid, Myeloid, Stromal, Immune)

merged_compartment_data <- left_join(sample.inf.final, wide_data_comp, by = "spatial_core_id")

long_data_compartment <- melt(
  merged_compartment_data,
  id.vars = c(
    "Study.ID", "spatial_core_id", "core_id", "stt_num", "tissue", "type",
    "subtype", "histological.subtype", "spatial.core.id.1", "Gender",
    "Grade.new", "Survival.status", "surv.time.new.CN.mod", "MCT.count",
    "MCT.status", "MCT.hilo", "neg", "MCT.low..1.9.", "tumor.size",
    "Mut.dect", "Gene.mut"
  ),
  measure.vars = c("Lymphoid", "Myeloid", "Stromal", "Immune"),
  variable.name = "compartment",
  value.name = "proportion"
)

# =========================================================
# 8. Pairwise comparison across all histological subtypes for compartments
# =========================================================
subtype_combinations_comp <- combn(unique(long_data_compartment$histological.subtype), 2, simplify = FALSE)

pairwise_subtype_comparison_compartments <- map_dfr(subtype_combinations_comp, function(pair) {
  s1 <- pair[1]
  s2 <- pair[2]

  long_data_compartment %>%
    filter(histological.subtype %in% c(s1, s2)) %>%
    group_by(compartment) %>%
    reframe({
      x <- proportion[histological.subtype == s1]
      y <- proportion[histological.subtype == s2]

      wt <- tryCatch(
        wilcox.test(x, y, exact = FALSE, correct = FALSE),
        error = function(e) NULL
      )

      W  <- if (is.null(wt)) NA_real_ else unname(wt$statistic)
      pv <- if (is.null(wt)) NA_real_ else wt$p.value

      cd <- suppressWarnings(tryCatch(
        as.numeric(effsize::cliff.delta(x, y, method = "asymptotic")$estimate),
        error = function(e) NA_real_
      ))

      fc <- (mean(x, na.rm = TRUE) + eps) / (mean(y, na.rm = TRUE) + eps)

      tibble(
        Fold_Change = fc,
        PValue = pv,
        W = W,
        CliffDelta = cd,
        group1 = s1,
        group2 = s2
      )
    }) %>%
    mutate(Comparison = paste0(group2, "_vs_", group1)) %>%
    select(compartment, Fold_Change, PValue, W, CliffDelta, Comparison)
}) %>%
  mutate(fdr = p.adjust(PValue, method = "fdr"))

write.csv(
  pairwise_subtype_comparison_compartments,
  file = file.path(pair_dir, "pairwise_subtype_comparison_compartments.csv"),
  row.names = FALSE
)

# =========================================================
# 9. LUAD vs LUSC comparison for compartments
# =========================================================
luad_lusc_compartment <- long_data_compartment %>%
  filter(subtype %in% c(cmp_g1, cmp_g2))

LUAD_vs_LUSC_compartment_comparison <- luad_lusc_compartment %>%
  group_by(compartment) %>%
  summarise(
    mean_g1 = mean(proportion[subtype == cmp_g1], na.rm = TRUE),
    mean_g2 = mean(proportion[subtype == cmp_g2], na.rm = TRUE),
    W = {
      x <- proportion[subtype == cmp_g2]
      y <- proportion[subtype == cmp_g1]
      tryCatch(unname(wilcox.test(x, y, exact = FALSE, correct = FALSE)$statistic),
               error = function(e) NA_real_)
    },
    PValue = {
      x <- proportion[subtype == cmp_g2]
      y <- proportion[subtype == cmp_g1]
      tryCatch(wilcox.test(x, y, exact = FALSE, correct = FALSE)$p.value,
               error = function(e) NA_real_)
    },
    CliffDelta = {
      x <- proportion[subtype == cmp_g2]
      y <- proportion[subtype == cmp_g1]
      suppressWarnings(tryCatch(
        as.numeric(effsize::cliff.delta(x, y, method = "asymptotic")$estimate),
        error = function(e) NA_real_
      ))
    },
    Fold_Change = (mean_g2 + eps) / (mean_g1 + eps),
    .groups = "drop"
  ) %>%
  mutate(
    fdr = p.adjust(PValue, method = "fdr"),
    Comparison = paste0(cmp_g1, "_vs_", cmp_g2)
  ) %>%
  transmute(
    Compartment = factor(compartment),
    Fold_Change, PValue, W, CliffDelta, Comparison, fdr
  )

write.csv(
  LUAD_vs_LUSC_compartment_comparison,
  file = file.path(pair_dir, "LUAD_vs_LUSC_compartment_comparison.csv"),
  row.names = FALSE
)

# =========================================================
# 10. Load clinical pairwise detail results and standardize outputs
# =========================================================
ade.pair.diff.results_df <- read.csv(
  file = file.path(detail_dir, "ade.pair.diff.results_df.csv"),
  header = TRUE, as.is = TRUE
)

squa.pair.diff.results_df <- read.csv(
  file = file.path(detail_dir, "squa.pair.diff.results_df.csv"),
  header = TRUE, as.is = TRUE
)

ade.pair.diff.results_df <- ade.pair.diff.results_df[, c("Cell_Type", "Fold_Change", "PValue", "W", "CliffDelta", "Comparison", "fdr")]
colnames(ade.pair.diff.results_df)[colnames(ade.pair.diff.results_df) == "Cell_Type"] <- "Term"
ade.pair.diff.results_df$Group <- "LUAD: Cell (clinical)"

squa.pair.diff.results_df$fdr <- p.adjust(squa.pair.diff.results_df$PValue, method = "fdr")
squa.pair.diff.results_df <- squa.pair.diff.results_df[, c("Cell_Type", "Fold_Change", "PValue", "W", "CliffDelta", "Comparison", "fdr")]
colnames(squa.pair.diff.results_df)[colnames(squa.pair.diff.results_df) == "Cell_Type"] <- "Term"
squa.pair.diff.results_df$Group <- "LUSC: Cell (clinical)"

pairwise_subtype_comparison_all_cells_std <- pairwise_subtype_comparison_all_cells[, c("Cell_Type", "Fold_Change", "PValue", "W", "CliffDelta", "Comparison", "fdr")]
colnames(pairwise_subtype_comparison_all_cells_std)[colnames(pairwise_subtype_comparison_all_cells_std) == "Cell_Type"] <- "Term"
pairwise_subtype_comparison_all_cells_std$Group <- "Cell (clinical)"

pairwise_subtype_comparison_compartments_std <- pairwise_subtype_comparison_compartments[, c("compartment", "Fold_Change", "PValue", "W", "CliffDelta", "Comparison", "fdr")]
colnames(pairwise_subtype_comparison_compartments_std)[colnames(pairwise_subtype_comparison_compartments_std) == "compartment"] <- "Term"
pairwise_subtype_comparison_compartments_std$Group <- "Compartment (clinical)"

LUAD_vs_LUSC_celltype_comparison_std <- LUAD_vs_LUSC_celltype_comparison[, c("Cell_Type", "Fold_Change", "PValue", "W", "CliffDelta", "Comparison", "fdr")]
colnames(LUAD_vs_LUSC_celltype_comparison_std)[colnames(LUAD_vs_LUSC_celltype_comparison_std) == "Cell_Type"] <- "Term"
LUAD_vs_LUSC_celltype_comparison_std$Group <- "Cell (LUAD vs LUSC)"

LUAD_vs_LUSC_compartment_comparison_std <- LUAD_vs_LUSC_compartment_comparison[, c("Compartment", "Fold_Change", "PValue", "W", "CliffDelta", "Comparison", "fdr")]
colnames(LUAD_vs_LUSC_compartment_comparison_std)[colnames(LUAD_vs_LUSC_compartment_comparison_std) == "Compartment"] <- "Term"
LUAD_vs_LUSC_compartment_comparison_std$Group <- "Compartment (LUAD vs LUSC)"

# =========================================================
# 11. Standardize comparison labels and drop Others
# =========================================================
rename_comparison_labels <- function(x) {
  x <- gsub("surv.time.new.CN.mod", "Surv_time", x, fixed = TRUE)
  x <- gsub("Grade.new", "Grade", x, fixed = TRUE)
  x <- gsub("MCT.status", "MCT", x, fixed = TRUE)
  x <- gsub("Mut.dect", "Mutation", x, fixed = TRUE)
  x <- gsub("tumor.size", "Tumor_size", x, fixed = TRUE)
  x
}

ade.pair.diff.results_df$Comparison                  <- rename_comparison_labels(ade.pair.diff.results_df$Comparison)
squa.pair.diff.results_df$Comparison                 <- rename_comparison_labels(squa.pair.diff.results_df$Comparison)
pairwise_subtype_comparison_all_cells_std$Comparison <- rename_comparison_labels(pairwise_subtype_comparison_all_cells_std$Comparison)
pairwise_subtype_comparison_compartments_std$Comparison <- rename_comparison_labels(pairwise_subtype_comparison_compartments_std$Comparison)
LUAD_vs_LUSC_celltype_comparison_std$Comparison      <- rename_comparison_labels(LUAD_vs_LUSC_celltype_comparison_std$Comparison)
LUAD_vs_LUSC_compartment_comparison_std$Comparison   <- rename_comparison_labels(LUAD_vs_LUSC_compartment_comparison_std$Comparison)

ade.pair.diff.results_df <- ade.pair.diff.results_df[ade.pair.diff.results_df$Term != "Others", ]
squa.pair.diff.results_df <- squa.pair.diff.results_df[squa.pair.diff.results_df$Term != "Others", ]
pairwise_subtype_comparison_all_cells_std <- pairwise_subtype_comparison_all_cells_std[pairwise_subtype_comparison_all_cells_std$Term != "Others", ]
LUAD_vs_LUSC_celltype_comparison_std <- LUAD_vs_LUSC_celltype_comparison_std[LUAD_vs_LUSC_celltype_comparison_std$Term != "Others", ]

# =========================================================
# 12. Combine final tables
# =========================================================
all.compare.cliffDelta.fdr.comb <- rbind.data.frame(
  ade.pair.diff.results_df,
  squa.pair.diff.results_df,
  pairwise_subtype_comparison_compartments_std,
  LUAD_vs_LUSC_celltype_comparison_std,
  LUAD_vs_LUSC_compartment_comparison_std
)

Term <- as.character(all.compare.cliffDelta.fdr.comb$Term)
Term[Term == "Endo"]        <- "EC"
Term[Term == "NK.like"]     <- "NK.cell"
Term[Term == "Plasma.cell"] <- "Plasma"
Term[Term == "Mono"]        <- "Monocyte"
Term[Term == "DCs"]         <- "DC"
all.compare.cliffDelta.fdr.comb$Term <- Term

write.csv(ade.pair.diff.results_df, file.path(final_dir, "ade.pair.diff.results_df.csv"), row.names = FALSE)
write.csv(squa.pair.diff.results_df, file.path(final_dir, "squa.pair.diff.results_df.csv"), row.names = FALSE)
write.csv(pairwise_subtype_comparison_all_cells_std, file.path(final_dir, "pairwise_subtype_comparison_all_cells.csv"), row.names = FALSE)
write.csv(pairwise_subtype_comparison_compartments_std, file.path(final_dir, "pairwise_subtype_comparison_compartments.csv"), row.names = FALSE)
write.csv(LUAD_vs_LUSC_celltype_comparison_std, file.path(final_dir, "LUAD_vs_LUSC_celltype_comparison.csv"), row.names = FALSE)
write.csv(LUAD_vs_LUSC_compartment_comparison_std, file.path(final_dir, "LUAD_vs_LUSC_compartment_comparison.csv"), row.names = FALSE)
write.csv(all.compare.cliffDelta.fdr.comb, file.path(final_dir, "all.compare.cliffDelta.fdr.comb.csv"), row.names = FALSE)

# =========================================================
# 13. Plot Cliff's delta with significance highlighting
# =========================================================
plot_delta_p <- function(df, title = "") {
  df1 <- df %>%
    mutate(
      sig = !is.na(PValue) & PValue < 0.05,
      sig_lab = factor(if_else(sig, "P < 0.05", "ns"), levels = c("ns", "P < 0.05")),
      ylab = if ("Comparison" %in% names(.)) paste0(Term, " | ", Comparison) else as.character(Term)
    ) %>%
    arrange(CliffDelta) %>%
    mutate(ylab = factor(ylab, levels = ylab))

  ggplot(df1, aes(x = CliffDelta, y = ylab, color = sig_lab)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_vline(
      xintercept = c(-0.147, 0.147, -0.33, 0.33, -0.474, 0.474),
      linetype = 3,
      linewidth = 0.3
    ) +
    geom_segment(
      aes(x = 0, xend = CliffDelta, y = ylab, yend = ylab),
      linewidth = 0.8,
      alpha = 0.6
    ) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("ns" = "grey40", "P < 0.05" = "#B2182B"), drop = FALSE) +
    coord_cartesian(xlim = c(-1, 1)) +
    labs(
      title = title,
      x = "Cliff's delta",
      y = NULL,
      color = NULL
    ) +
    theme_classic(base_size = 12)
}

p1 <- plot_delta_p(ade.pair.diff.results_df, "LUAD: Cell level (pairwise)")
p2 <- plot_delta_p(squa.pair.diff.results_df, "LUSC: Cell level (pairwise)")

combined_plot <- plot_grid(p1, p2, nrow = 1, align = "v")

ggsave(
  filename = file.path(plot_dir, "ade.squa.combined_bubble_plots.pdf"),
  plot = combined_plot,
  width = 18.41,
  height = 5.77
)