rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(purrr)
  library(ggplot2)
  library(grid)
})

# =========================================================
# 1. Define paths
# =========================================================
base_dir <- "/cell.anno"
dens_dir <- file.path(base_dir, "cell.density.res")
meta_dir <- file.path(base_dir, "/results/cell.fraction.diff.res")
out_dir  <- file.path(base_dir, "/results/TMA_cell.freq.corr.res")
plot_dir <- file.path(out_dir, "plot")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. Load Seurat object and sample metadata
# =========================================================
codex.obj <- readRDS(file.path(dens_dir, "codex.obj.update.337.rds"))
filtered_meta_data <- codex.obj@meta.data

sample.inf.comb <- read.csv(
  file.path(meta_dir, "sample.inf.comb.csv"),
  header = TRUE,
  as.is = TRUE
)
sample.inf.comb <- sample.inf.comb[, -1]

colnames(sample.inf.comb)[colnames(sample.inf.comb) == "Study.ID"] <- "core_id"
colnames(sample.inf.comb)[colnames(sample.inf.comb) == "subtype.y"] <- "subtype"
colnames(sample.inf.comb)[colnames(sample.inf.comb) == "spatial_core_id"] <- "spatial.core.id"

sample.inf <- sample.inf.comb

his.subtype.new <- sample.inf$histological.subtype
his.subtype.new[his.subtype.new %in% c("Lepidic", "Papillary", "Acinar")] <- "ade.well"
his.subtype.new[his.subtype.new %in% c("Micropapillary", "Solid")] <- "ade.poor"
his.subtype.new[his.subtype.new %in% c("Keratinizing squamous cell carcinoma")] <- "squa.poor"
his.subtype.new[his.subtype.new %in% c("Non-keratinizing squamous cell carcinoma")] <- "squa.well"
sample.inf$his.subtype.new <- his.subtype.new

sample.inf <- sample.inf %>%
  dplyr::select(spatial.core.id, subtype, his.subtype.new, Grade.new)

# =========================================================
# 3. Compute core-level cell-type abundance table
# =========================================================
meta.data <- filtered_meta_data

counts <- meta.data %>%
  group_by(spatial.core.id, cell.anno.3st) %>%
  summarise(count = n(), .groups = "drop")

totals <- counts %>%
  group_by(spatial.core.id) %>%
  summarise(total = sum(count), .groups = "drop")

proportions <- counts %>%
  left_join(totals, by = "spatial.core.id") %>%
  mutate(proportion = count / total)

proportions_adjusted <- dplyr::select(proportions, -proportion)

wide_data <- proportions_adjusted %>%
  pivot_wider(
    names_from = cell.anno.3st,
    values_from = count,
    values_fill = list(count = 0)
  ) %>%
  select(-total)

# =========================================================
# 4. Define cell-type order and colors
# =========================================================
color.value.anno <- c(
  "Cancer.cell"   = "#B17A86",
  "Tcyto"         = "#6BD089",
  "Th"            = "#A5CA6B",
  "Treg"          = "#D8E063",
  "NK.like"       = "#138535",
  "Neutrophil"    = "#00D5F2",
  "Plasma.cell"   = "#75C5D9",
  "Bcell"         = "#2DA2BF",
  "Mono"          = "#7A297B",
  "M1"            = "#C785C8",
  "M2"            = "#BDB6D6",
  "DCs"           = "#586B8C",
  "CAF"           = "#FFCB00",
  "Smooth.muscle" = "#E19F73",
  "Endo"          = "#F08B1E",
  "Others"        = "gray90"
)

cell_types <- names(color.value.anno)

# =========================================================
# 5. Helper function for pairwise Spearman correlation
# =========================================================
compute_correlation <- function(df, group_col) {
  combinations <- expand.grid(cell_types, cell_types, stringsAsFactors = FALSE)
  colnames(combinations) <- c("from", "to")

  corr_list <- lapply(seq_len(nrow(combinations)), function(i) {
    from_cell <- combinations$from[i]
    to_cell   <- combinations$to[i]

    test_result <- tryCatch(
      cor.test(df[[from_cell]], df[[to_cell]], method = "spearman"),
      error = function(e) NULL
    )

    data.frame(
      from = from_cell,
      to = to_cell,
      rho = if (is.null(test_result)) NA_real_ else unname(test_result$estimate),
      pvalue = if (is.null(test_result)) NA_real_ else test_result$p.value
    )
  })

  bind_rows(corr_list)
}

# =========================================================
# 6. Subtype-level co-occurrence analysis
# =========================================================
joined_data_subtype <- left_join(sample.inf, wide_data, by = "spatial.core.id")

wide_frac_subtype <- joined_data_subtype %>%
  dplyr::select(spatial.core.id, subtype, all_of(cell_types))

cor_results_subtype <- wide_frac_subtype %>%
  group_by(subtype) %>%
  nest() %>%
  mutate(cor_data = map(data, ~ compute_correlation(.x, "subtype"))) %>%
  unnest(cor_data) %>%
  ungroup()

cor_results_subtype <- cor_results_subtype %>%
  group_by(subtype) %>%
  mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()

write.csv(
  cor_results_subtype,
  file = file.path(out_dir, "cell_fraction_correlation_by_subtype.csv"),
  row.names = FALSE
)

# Prepare plotting data
cor_results_subtype$from <- factor(cor_results_subtype$from, levels = cell_types)
cor_results_subtype$to <- factor(cor_results_subtype$to, levels = rev(cell_types))
cor_results_subtype$subtype <- factor(cor_results_subtype$subtype, levels = rev(unique(cor_results_subtype$subtype)))

cor_results_subtype <- cor_results_subtype %>%
  mutate(
    y_pos = as.numeric(factor(to)),
    subtype_offset = ifelse(subtype == "adeno", 0.5, 0),
    label_pos = y_pos + 0.25
  )

p_subtype <- ggplot(cor_results_subtype, aes(x = from, y = y_pos + subtype_offset)) +
  geom_tile(aes(fill = rho), height = 0.5, width = 0.9) +
  geom_tile(
    aes(y = label_pos),
    height = 1,
    width = 1,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    space = "Lab"
  ) +
  scale_y_continuous(
    breaks = unique(cor_results_subtype$label_pos),
    labels = unique(cor_results_subtype$to)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.spacing = unit(0, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.y = element_blank()
  )

pdf(file.path(plot_dir, "TMA_cell.frac_subtype_corr_heatmap.subtype.pdf"), width = 7.52, height = 6.47)
print(p_subtype)
dev.off()

# =========================================================
# 7. Histologic-group co-occurrence analysis
# =========================================================
joined_data_his <- left_join(sample.inf, wide_data, by = "spatial.core.id")

wide_frac_his <- joined_data_his %>%
  dplyr::select(spatial.core.id, his.subtype.new, all_of(cell_types))

cor_results_his <- wide_frac_his %>%
  group_by(his.subtype.new) %>%
  nest() %>%
  mutate(cor_data = map(data, ~ compute_correlation(.x, "his.subtype.new"))) %>%
  unnest(cor_data) %>%
  ungroup()

cor_results_his <- cor_results_his %>%
  group_by(his.subtype.new) %>%
  mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()

write.csv(
  cor_results_his,
  file = file.path(out_dir, "cell_fraction_correlation_by_histologic_group.csv"),
  row.names = FALSE
)

levels_to_set <- rev(c("ade.well", "ade.poor", "squa.well", "squa.poor"))
cor_results_his$his.subtype.new <- factor(cor_results_his$his.subtype.new, levels = levels_to_set)

cor_results_his <- cor_results_his %>%
  mutate(
    from = factor(from, levels = cell_types),
    to = factor(to, levels = rev(cell_types)),
    y_pos = as.numeric(to),
    x_pos = as.numeric(from),
    offset = as.numeric(factor(his.subtype.new, levels = levels_to_set)) / 4
  )

p_his <- ggplot(cor_results_his, aes(x = x_pos, y = y_pos + offset)) +
  geom_tile(aes(fill = rho), height = 1 / 4, width = 1, color = NA) +
  geom_tile(
    aes(y = y_pos + 0.5 + 1 / 12),
    fill = NA,
    color = "black",
    linewidth = 0.5,
    height = 1,
    width = 1
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    space = "Lab",
    name = "Spearman rho"
  ) +
  scale_y_continuous(
    breaks = unique(cor_results_his$y_pos + 0.5),
    labels = unique(cor_results_his$to)
  ) +
  scale_x_continuous(
    breaks = unique(cor_results_his$x_pos),
    labels = unique(cor_results_his$from)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.spacing = unit(0, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.y = element_blank()
  )

pdf(file.path(plot_dir, "TMA_cell.frac_his.subtype_corr_heatmap.subtype.pdf"), width = 7.52, height = 6.47)
print(p_his)
dev.off()