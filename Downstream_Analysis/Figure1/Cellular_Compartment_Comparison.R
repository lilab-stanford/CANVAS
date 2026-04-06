rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(data.table)
  library(Seurat)
  library(ggbeeswarm)
})

# =========================================================
# 1. Define paths
# =========================================================
base_dir <- "/cell.anno/cell.density.res"
subtype_dir <- file.path(base_dir, "histological.subtype")
out_dir <- file.path(base_dir, "path.subtype.diff/results/cell.fraction.diff.res")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

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
  select(spatial.core.id, core_id, study_id, stt_num, tissue, type, subtype,
         spatial.core.id.1, each.core.id) %>%
  distinct(study_id, .keep_all = TRUE)

rownames(meta.unique) <- meta.unique$spatial.core.id
colnames(meta.unique)[colnames(meta.unique) == "study_id"] <- "Study.ID"

sample.unique <- sample.inf %>%
  filter(!is.na(tissue)) %>%
  select(spatial.core.id, core_id, Study.ID, stt_num, tissue, type, subtype,
         spatial.core.id.1, each.core.id, Age, Gender, Survival.status,
         surv.time.new.CN.mod, histological.subtype, Grade.new, Stage,
         Invasion, MCT.status, MCT.hilo, tumor.size, Mut.dect, Gene.mut) %>%
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

fwrite(merged_data, file = file.path(out_dir, "merged_data.csv"))

# =========================================================
# 5. Convert to long format for subtype comparison
# =========================================================
cell_types <- c(
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
  measure.vars = intersect(cell_types, colnames(merged_data)),
  variable.name = "cell_type",
  value.name = "proportion"
)

long_data <- long_data %>%
  filter(subtype %in% c("adeno", "squamous cell"))

# =========================================================
# 6. Compute summary statistics and Wilcoxon tests
# =========================================================
summary_stats <- long_data %>%
  group_by(cell_type, subtype) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop")

summary_wide <- summary_stats %>%
  pivot_wider(names_from = subtype, values_from = mean_prop) %>%
  mutate(fold_change = adeno / `squamous cell`)

pvals <- long_data %>%
  group_by(cell_type) %>%
  summarise(
    p_value = wilcox.test(proportion ~ subtype, exact = FALSE)$p.value,
    .groups = "drop"
  )

plot_data <- left_join(summary_wide, pvals, by = "cell_type") %>%
  arrange(p_value)

fwrite(plot_data, file = file.path(out_dir, "plot_data.csv"))

# =========================================================
# 7. Generate the main subtype comparison plot
# =========================================================
cell_order <- c(
  "Others", "NK.like", "DCs", "Smooth.muscle", "Plasma.cell", "Bcell",
  "Treg", "Mono", "M1", "Neutrophil", "M2", "Tcyto", "CAF", "Endo",
  "Th", "Cancer.cell"
)

color.value.anno <- c(
  "Cancer.cell"  = "#B17A86",
  "Tcyto"        = "#6BD089",
  "Th"           = "#A5CA6B",
  "Treg"         = "#D8E063",
  "NK.like"      = "#138535",
  "Neutrophil"   = "#00D5F2",
  "CAF"          = "#FFCB00",
  "Smooth.muscle"= "#E19F73",
  "Endo"         = "#F08B1E",
  "Plasma.cell"  = "#75C5D9",
  "Bcell"        = "#2DA2BF",
  "Mono"         = "#7A297B",
  "M1"           = "#C785C8",
  "M2"           = "#BDB6D6",
  "DCs"          = "#586B8C",
  "Others"       = "gray90"
)

plot_df <- long_data %>%
  mutate(
    cell_type = factor(cell_type, levels = cell_order),
    subtype   = factor(subtype, levels = c("adeno", "squamous cell")),
    log_proportion = log2(proportion + 1)
  )

BAR_WIDTH <- 0.42
BAR_SIZE  <- 0.7
DOT_SIZE  <- 0.2
DOT_ALPHA <- 0.6
JIT_BW    <- 0.08
JIT_DW    <- 0.10

p_all <- ggplot(plot_df, aes(x = cell_type, y = log_proportion, color = cell_type)) +
  geom_quasirandom(
    dodge.width = JIT_DW,
    bandwidth = JIT_BW,
    size = DOT_SIZE,
    alpha = DOT_ALPHA,
    varwidth = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "errorbar",
    aes(ymin = after_stat(y), ymax = after_stat(y)),
    width = BAR_WIDTH,
    linewidth = BAR_SIZE,
    color = "black"
  ) +
  coord_flip() +
  facet_wrap(~ subtype, nrow = 1, scales = "free_y") +
  scale_color_manual(values = color.value.anno, drop = FALSE) +
  labs(
    x = "Cell type",
    y = "Proportion of total cells [log2(% + 1)]",
    title = "Cell-type composition by subtype"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

ggsave(
  filename = file.path(out_dir, "cell_type_proportion_by_subtype_log2.pdf"),
  plot = p_all,
  width = 5.81,
  height = 4.73
)

# =========================================================
# 8. Save core objects and processed tables
# =========================================================
saveRDS(TMA.codex.obj, file = file.path(out_dir, "TMA.codex.obj.rds"))
fwrite(meta.data, file = file.path(out_dir, "filtered_meta_data.csv"))
fwrite(long_data, file = file.path(out_dir, "long_data.csv"))