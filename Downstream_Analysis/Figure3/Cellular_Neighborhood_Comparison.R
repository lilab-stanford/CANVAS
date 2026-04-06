rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

# =========================================================
# 1. Define paths
# =========================================================
base_dir <- "/cell.anno"
cn_file  <- file.path(base_dir, "cell.density.res/sample.inf.csv")
sam_file <- file.path(base_dir, "cell.density.res/histological.subtype/sample.inf.comb.csv")
mut_file <- file.path(base_dir, "cell.density.res/histological.subtype/tma_mut.csv")

result_dir <- file.path(base_dir, "/CN.sum.results_new/results")
plot_dir   <- file.path(base_dir, "/CN.sum.results_new/plot")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

eps <- 1e-12

# =========================================================
# 2. Load CN proportions and sample metadata
# =========================================================
cn_data <- read.csv(cn_file, header = TRUE, row.names = 1, as.is = TRUE)
wide_proportions <- cn_data[, c("spatial.core.id", paste0("CN", sprintf("%02d", 1:10)))]

sample.inf <- read.csv(sam_file, header = TRUE, row.names = 1, as.is = TRUE)
sample.inf <- sample.inf[!is.na(sample.inf$subtype), ]
sample.inf[sample.inf$subtype == "squamous cell", "histological.subtype"] <-
  sample.inf[sample.inf$subtype == "squamous cell", "histological.subtype.pre"]

drop_cols <- intersect(colnames(sample.inf), c("EGFR", "ALK", "ROS1", "KRAS"))
if (length(drop_cols) > 0) {
  sample.inf <- sample.inf[, !colnames(sample.inf) %in% drop_cols, drop = FALSE]
}

mut.info <- read.csv(mut_file, header = TRUE, row.names = 1, as.is = TRUE)
mut.info <- mut.info[, c("MRN", "Donor.ID", "tumor.size", "Mut.dect", "Gene.mut")]
sample.inf <- merge(sample.inf, mut.info, by.x = "Donor.no.", by.y = "Donor.ID", all = TRUE)

# =========================================================
# 3. Prepare clinical variables for comparison
# =========================================================
prepare_clinical_data <- function(df, subtype_name) {
  dat <- df[df$subtype == subtype_name, , drop = FALSE]

  dat <- dat %>%
    mutate(
      Age = ifelse(Age > 70, "high", "low"),
      Gender = ifelse(Gender == "male", "high", "low"),
      surv.time = surv.time.new.CN.mod,
      surv.time.new.CN.mod = ifelse(
        surv.time.new.CN.mod > median(surv.time.new.CN.mod, na.rm = TRUE),
        "high", "low"
      ),
      Grade.new = ifelse(Grade.new %in% c("well", "moderate"), "high", "low"),
      Stage = ifelse(Stage %in% c("Stage3", "Stage4"), "high", "low"),
      Invasion = ifelse(Invasion == "", "low", "high"),
      MCT.status = ifelse(MCT.status %in% c("high", "intermed"), "high", "low"),
      tumor.size = ifelse(tumor.size > 3, "high", "low"),
      Mut.dect = ifelse(Mut.dect == "TRUE", "low", NA)
    )

  dat$Mut.dect[which(dat$Gene.mut == "TRUE")] <- "high"
  dat
}

add_histology_pairwise_columns <- function(dat, subtype_levels) {
  dat$histological.subtype <- factor(dat$histological.subtype, levels = subtype_levels)
  combinations <- combn(levels(dat$histological.subtype), 2, simplify = FALSE)

  for (combo in combinations) {
    col_name <- paste0(combo[1], "_vs_", combo[2])
    dat[[col_name]] <- ifelse(
      dat$histological.subtype == combo[1], "low",
      ifelse(dat$histological.subtype == combo[2], "high", NA)
    )
  }
  dat
}

# =========================================================
# 4. Run CN vs parameter comparison
# =========================================================
run_cn_parameter_test <- function(meta_df, cn_df, comparison_columns) {
  inter_ids <- intersect(meta_df$spatial.core.id, cn_df$spatial.core.id)
  cn_sub <- cn_df[cn_df$spatial.core.id %in% inter_ids, , drop = FALSE]
  merged_data <- merge(meta_df, cn_sub, by = "spatial.core.id", all = TRUE)

  cn_cols <- setdiff(colnames(cn_sub), "spatial.core.id")
  results_df <- data.frame()

  for (comp in comparison_columns) {
    comp_results <- data.frame(Cell_Type = cn_cols, Fold_Change = NA, PValue = NA)

    for (cn in cn_cols) {
      high_data <- merged_data %>% filter(.data[[comp]] == "high") %>% pull(.data[[cn]])
      low_data  <- merged_data %>% filter(.data[[comp]] == "low")  %>% pull(.data[[cn]])

      if (length(high_data) > 0 && length(low_data) > 0) {
        wt <- tryCatch(wilcox.test(high_data, low_data, exact = FALSE), error = function(e) NULL)

        if (!is.null(wt)) {
          fc <- (mean(high_data, na.rm = TRUE) + eps) / (mean(low_data, na.rm = TRUE) + eps)
          comp_results[comp_results$Cell_Type == cn, c("Fold_Change", "PValue")] <- c(fc, wt$p.value)
        }
      }
    }

    comp_results$Comparison <- comp
    results_df <- rbind(results_df, comp_results)
  }

  results_df
}

# =========================================================
# 5. Format bubble plot data
# =========================================================
prepare_bubble_data <- function(results_df, comparison_levels) {
  results_df$logp <- -log10(pmax(results_df$PValue, .Machine$double.xmin))

  results_df$Size <- ifelse(
    results_df$PValue <= 0.001, 3,
    ifelse(
      results_df$PValue <= 0.01, 2,
      ifelse(
        results_df$PValue <= 0.05, 1,
        pmin(results_df$logp / max(results_df$logp, na.rm = TRUE), 1)
      )
    )
  )

  results_df$para.group.anno <- ifelse(
    results_df$Fold_Change > 1,
    paste0(results_df$Comparison, "_high"),
    paste0(results_df$Comparison, "_low")
  )

  cn_levels <- paste0("CN", sprintf("%02d", 1:10))
  results_df$Cell_Type <- factor(results_df$Cell_Type, levels = cn_levels)
  results_df$Comparison <- factor(results_df$Comparison, levels = rev(comparison_levels))
  results_df
}

plot_cn_bubble <- function(results_df, fill_colors, file_out, width, height) {
  p <- ggplot(results_df, aes(x = Cell_Type, y = Comparison, size = Size, fill = para.group.anno)) +
    geom_point(shape = 21, colour = "black") +
    scale_size_continuous(
      breaks = c(3, 2, 1),
      labels = c("P <= 0.001", "P <= 0.01", "P <= 0.05"),
      range = c(1, 12)
    ) +
    scale_fill_manual(values = fill_colors, breaks = names(fill_colors), drop = FALSE) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(),
      legend.position = "right"
    ) +
    labs(
      x = "CN",
      y = "Comparison",
      title = "Bubble Plot of CN Comparisons",
      subtitle = "Bubble size indicates P-value significance; fill color indicates direction of fold change"
    ) +
    guides(
      size = guide_legend(
        order = 1,
        title = "P-value Significance",
        byrow = TRUE,
        nrow = 3,
        ncol = 1,
        override.aes = list(fill = "black")
      ),
      fill = guide_legend(
        order = 2,
        title = "Parameter Group",
        byrow = TRUE,
        ncol = 2
      )
    )

  ggsave(file_out, p, width = width, height = height)
  p
}

# =========================================================
# 6. Define color palette
# =========================================================
detail_fill_colors <- c(
  "Age_low" = "#9090D5",
  "Age_high" = "#00008B",
  "Gender_low" = "#D7A0BC",
  "Gender_high" = "#FF69B4",
  "surv.time.new.CN.mod_low" = "#8CD58C",
  "surv.time.new.CN.mod_high" = "#008000",
  "Grade.new_low" = "#D29984",
  "Grade.new_high" = "#FF4500",
  "Stage_low" = "#CE9191",
  "Stage_high" = "#A52A2A",
  "Invasion_low" = "#CDA88E",
  "Invasion_high" = "#8B4513",
  "MCT.status_low" = "#D8D6B4",
  "MCT.status_high" = "#BDB76B",
  "tumor.size_low" = "#B19BCD",
  "tumor.size_high" = "#553879",
  "Mut.dect_low" = "#CFA085",
  "Mut.dect_high" = "#B0480F",
  "Lepidic_vs_Papillary_low" = "#EDD3ED",
  "Lepidic_vs_Papillary_high" = "#DAAADB",
  "Lepidic_vs_Acinar_low" = "#EDD3ED",
  "Lepidic_vs_Acinar_high" = "#DAAADB",
  "Lepidic_vs_Micropapillary_low" = "#EDD3ED",
  "Lepidic_vs_Micropapillary_high" = "#AA94D7",
  "Lepidic_vs_Solid_low" = "#EDD3ED",
  "Lepidic_vs_Solid_high" = "#9932CC",
  "Papillary_vs_Solid_low" = "#DAAADB",
  "Papillary_vs_Solid_high" = "#9932CC",
  "Acinar_vs_Solid_low" = "#DAAADB",
  "Acinar_vs_Solid_high" = "#9932CC",
  "Micropapillary_vs_Solid_low" = "#AA94D7",
  "Micropapillary_vs_Solid_high" = "#9932CC",
  "Non-kerat_vs_Kerat_low" = "#EDD3ED",
  "Non-kerat_vs_Kerat_high" = "#BD84BA"
)

# =========================================================
# 7. LUAD analysis
# =========================================================
luad_meta <- prepare_clinical_data(sample.inf, "adeno")
luad_meta <- add_histology_pairwise_columns(
  luad_meta,
  subtype_levels = c("Lepidic", "Papillary", "Acinar", "Micropapillary", "Solid")
)

luad_comparisons <- c(
  "Age", "Gender", "surv.time.new.CN.mod", "Grade.new", "Stage", "Invasion",
  "MCT.status", "tumor.size", "Mut.dect",
  "Lepidic_vs_Papillary", "Lepidic_vs_Acinar", "Lepidic_vs_Micropapillary",
  "Lepidic_vs_Solid", "Papillary_vs_Acinar", "Papillary_vs_Micropapillary",
  "Papillary_vs_Solid", "Acinar_vs_Micropapillary", "Acinar_vs_Solid",
  "Micropapillary_vs_Solid"
)

luad_results <- run_cn_parameter_test(luad_meta, wide_proportions, luad_comparisons)

drop_luad <- c("Papillary_vs_Acinar", "Papillary_vs_Micropapillary", "Acinar_vs_Micropapillary")
luad_results <- luad_results[!luad_results$Comparison %in% drop_luad, ]

write.csv(luad_results, file.path(result_dir, "ade.pair.diff.results_df.csv"), row.names = FALSE)

luad_plot_levels <- c(
  "Age", "Gender", "surv.time.new.CN.mod", "Grade.new", "Stage", "Invasion",
  "MCT.status", "tumor.size", "Mut.dect",
  "Lepidic_vs_Papillary", "Lepidic_vs_Acinar", "Lepidic_vs_Micropapillary",
  "Lepidic_vs_Solid", "Papillary_vs_Solid", "Acinar_vs_Solid",
  "Micropapillary_vs_Solid"
)

luad_results_plot <- prepare_bubble_data(luad_results, luad_plot_levels)

ade.bubble.plot <- plot_cn_bubble(
  luad_results_plot,
  fill_colors = detail_fill_colors,
  file_out = file.path(plot_dir, "adeno_CN_bubble_plots.pdf"),
  width = 11.5,
  height = 7.0
)

# =========================================================
# 8. LUSC analysis
# =========================================================
lusc_meta <- prepare_clinical_data(sample.inf, "squamous cell")
lusc_meta$histological.subtype[lusc_meta$histological.subtype == "Non-keratinizing squamous cell carcinoma"] <- "Non-kerat"
lusc_meta$histological.subtype[lusc_meta$histological.subtype == "Keratinizing squamous cell carcinoma"] <- "Kerat"

lusc_meta <- add_histology_pairwise_columns(
  lusc_meta,
  subtype_levels = c("Non-kerat", "Kerat")
)

lusc_comparisons <- c(
  "Age", "Gender", "surv.time.new.CN.mod", "Grade.new", "Stage",
  "Invasion", "MCT.status", "tumor.size", "Mut.dect", "Non-kerat_vs_Kerat"
)

lusc_results <- run_cn_parameter_test(lusc_meta, wide_proportions, lusc_comparisons)

write.csv(lusc_results, file.path(result_dir, "squa.pair.diff.results_df.csv"), row.names = FALSE)

lusc_plot_levels <- c(
  "Age", "Gender", "surv.time.new.CN.mod", "Grade.new", "Stage",
  "Invasion", "MCT.status", "tumor.size", "Mut.dect", "Non-kerat_vs_Kerat"
)

lusc_results_plot <- prepare_bubble_data(lusc_results, lusc_plot_levels)

squa.bubble.plot <- plot_cn_bubble(
  lusc_results_plot,
  fill_colors = detail_fill_colors,
  file_out = file.path(plot_dir, "squa_CN_bubble_plots.pdf"),
  width = 11.0,
  height = 6.8
)

# =========================================================
# 9. Combined panel
# =========================================================
combined_plot <- plot_grid(ade.bubble.plot, squa.bubble.plot, nrow = 1, align = "v")

ggsave(
  filename = file.path(plot_dir, "ade.squa.combined_CN_bubble_plots.pdf"),
  plot = combined_plot,
  width = 22,
  height = 6.2
)