rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# =========================================================
# 1. Define paths and create output directories
# =========================================================
base_dir <- "/Pra_fdr_bubble_plot"
result_dir <- "/Para_bubble_plot/results"
plot_dir <- file.path(base_dir, "plot")

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "RCode"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "rawdata"), showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 2. Define color palettes
# =========================================================
cell_type_colors <- c(
  "Cancer.cell"   = "#B17A86",
  "Tcyto"         = "#6BD089",
  "Th"            = "#A5CA6B",
  "Treg"          = "#D8E063",
  "NK.like"       = "#138535",
  "Neutrophil"    = "#00D5F2",
  "CAF"           = "#FFCB00",
  "Smooth.muscle" = "#E19F73",
  "Endo"          = "#F08B1E",
  "Plasma.cell"   = "#75C5D9",
  "Bcell"         = "#2DA2BF",
  "Mono"          = "#7A297B",
  "M1"            = "#C785C8",
  "M2"            = "#BDB6D6",
  "DCs"           = "#586B8C",
  "Others"        = "gray90"
)

detail_fill_colors <- c(
  "Age_low"                     = "#9090D5",
  "Age_high"                    = "#00008B",
  "Gender_low"                  = "#D7A0BC",
  "Gender_high"                 = "#FF69B4",
  "surv.time.new.CN.mod_low"    = "#8CD58C",
  "surv.time.new.CN.mod_high"   = "#008000",
  "Grade.new_low"               = "#D29984",
  "Grade.new_high"              = "#FF4500",
  "Stage_low"                   = "#CE9191",
  "Stage_high"                  = "#A52A2A",
  "Invasion_low"                = "#CDA88E",
  "Invasion_high"               = "#8B4513",
  "MCT.status_low"              = "#D8D6B4",
  "MCT.status_high"             = "#BDB76B",
  "tumor.size_low"              = "#B19BCD",
  "tumor.size_high"             = "#553879",
  "Mut.dect_low"                = "#CFA085",
  "Mut.dect_high"               = "#B0480F",
  "Lepidic_vs_Papillary_low"    = "#EDD3ED",
  "Lepidic_vs_Papillary_high"   = "#DAAADB",
  "Lepidic_vs_Acinar_low"       = "#EDD3ED",
  "Lepidic_vs_Acinar_high"      = "#DAAADB",
  "Lepidic_vs_Micropapillary_low"  = "#EDD3ED",
  "Lepidic_vs_Micropapillary_high" = "#AA94D7",
  "Lepidic_vs_Solid_low"        = "#EDD3ED",
  "Lepidic_vs_Solid_high"       = "#9932CC",
  "Papillary_vs_Solid_low"      = "#DAAADB",
  "Papillary_vs_Solid_high"     = "#9932CC",
  "Acinar_vs_Solid_low"         = "#DAAADB",
  "Acinar_vs_Solid_high"        = "#9932CC",
  "Micropapillary_vs_Solid_low" = "#AA94D7",
  "Micropapillary_vs_Solid_high"= "#9932CC",
  "Non-kerat_vs_Kerat_low"      = "#FABF8F",
  "Non-kerat_vs_Kerat_high"     = "#B96927"
)

# =========================================================
# 3. Define shared helper functions
# =========================================================
prepare_results_df <- function(df) {
  df <- df %>%
    group_by(Comparison) %>%
    mutate(
      FDR = p.adjust(PValue, method = "BH")
    ) %>%
    ungroup()

  df$PValue <- df$FDR
  df$logp <- -log10(pmax(df$PValue, .Machine$double.xmin))

  df$Size <- case_when(
    df$PValue <= 0.01 ~ 2.5,
    df$PValue <= 0.10 ~ 1.5,
    df$PValue <= 0.20 ~ 0.8,
    df$PValue <= 0.30 ~ 0.6,
    df$PValue <= 0.40 ~ 0.5,
    df$PValue <= 0.70 ~ 0.4,
    df$PValue <= 1.00 ~ 0.3,
    TRUE ~ pmax(df$logp / max(df$logp, na.rm = TRUE), 0.3)
  )

  df$para.group.anno <- ifelse(
    df$Fold_Change > 1,
    paste0(df$Comparison, "_high"),
    paste0(df$Comparison, "_low")
  )

  df$Cell_Type <- factor(df$Cell_Type, levels = names(cell_type_colors))
  df
}

make_bubble_plot <- function(results_df, comparison_levels, width = 11, height = 7) {
  results_df$Comparison <- factor(results_df$Comparison, levels = rev(comparison_levels))

  p <- ggplot(
    results_df,
    aes(x = Cell_Type, y = Comparison, size = Size, fill = para.group.anno)
  ) +
    geom_point(shape = 21, colour = "black") +
    scale_size_continuous(
      breaks = c(2.5, 1.5, 0.8, 0.6, 0.5, 0.4, 0.3),
      labels = c("P ≤ 0.01", "P ≤ 0.10", "P ≤ 0.20", "P ≤ 0.30", "P ≤ 0.40", "P ≤ 0.70", "P ≤ 1.00"),
      range = c(1.5, 10)
    ) +
    scale_fill_manual(values = detail_fill_colors, breaks = names(detail_fill_colors), drop = FALSE) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(),
      legend.position = "right"
    ) +
    labs(
      x = "Comparison",
      y = "Cell Type",
      title = "Bubble Plot of Cell Type Comparisons",
      subtitle = "Bubble size indicates adjusted significance and fill color indicates direction of fold change"
    ) +
    guides(
      size = guide_legend(
        order = 1,
        title = "Adjusted P-value",
        byrow = TRUE,
        override.aes = list(fill = "black")
      ),
      fill = guide_legend(
        order = 2,
        title = "Parameter Group",
        byrow = TRUE,
        ncol = 2
      )
    )

  return(p)
}

# =========================================================
# 4. LUAD plot
# =========================================================
ade_results <- read.csv(
  file = file.path(result_dir, "ade.pair.diff.results_df.csv"),
  header = TRUE,
  as.is = TRUE,
  row.names = 1
)

ade_results <- prepare_results_df(ade_results)

ade_drop <- c("Papillary_vs_Acinar", "Papillary_vs_Micropapillary", "Acinar_vs_Micropapillary")
ade_results <- ade_results[!ade_results$Comparison %in% ade_drop, ]

ade_comparison_levels <- c(
  "Age", "Gender", "surv.time.new.CN.mod", "Grade.new", "Stage",
  "Invasion", "MCT.status", "tumor.size", "Mut.dect",
  "Lepidic_vs_Papillary", "Lepidic_vs_Acinar", "Lepidic_vs_Micropapillary",
  "Lepidic_vs_Solid", "Papillary_vs_Solid", "Acinar_vs_Solid",
  "Micropapillary_vs_Solid"
)

ade_bubble_plot <- make_bubble_plot(
  results_df = ade_results,
  comparison_levels = ade_comparison_levels
)

pdf(file = file.path(plot_dir, "adeno_cell.frac_bubble_Plots.pdf"), width = 11.92, height = 7.47)
print(ade_bubble_plot)
dev.off()

# =========================================================
# 5. LUSC plot
# =========================================================
squa_results <- read.csv(
  file = file.path(result_dir, "squa.pair.diff.results_df.csv"),
  header = TRUE,
  as.is = TRUE,
  row.names = 1
)

squa_results <- prepare_results_df(squa_results)

squa_comparison_levels <- c(
  "Age", "Gender", "surv.time.new.CN.mod", "Grade.new", "Stage",
  "Invasion", "MCT.status", "tumor.size", "Mut.dect",
  "Non-kerat_vs_Kerat"
)

squa_bubble_plot <- make_bubble_plot(
  results_df = squa_results,
  comparison_levels = squa_comparison_levels
)

pdf(file = file.path(plot_dir, "squ_cell.frac_bubble_Plots.pdf"), width = 11.07, height = 7.28)
print(squa_bubble_plot)
dev.off()

# =========================================================
# 6. Combined plot
# =========================================================
combined_plot <- plot_grid(ade_bubble_plot, squa_bubble_plot, nrow = 1, align = "v")

ggsave(
  filename = file.path(plot_dir, "ade.squa.combined_bubble_plots.pdf"),
  plot = combined_plot,
  width = 22,
  height = 6.5
)