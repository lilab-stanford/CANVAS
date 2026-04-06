#!/usr/bin/env Rscript
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pROC)
  library(ggsci)
  library(patchwork)
})

# =========================================================
# 1. Define output directory
# =========================================================
out_dir <- "/model_benchmarked/results_new_compare/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 2. Load prediction results
# =========================================================
# CANVAS
val_canvas_res   <- fread("/val.csv")
test_canvas_res  <- fread("/tma2.csv")

# uni
val_uni_res      <- fread("/val_uni.csv")
test_uni_res     <- fread("/tma2_uni.csv")

# virchow
val_virchow_res  <- fread("/val_virchow.csv")
test_virchow_res <- fread("/tma2_virchow.csv")

# resnet
val_resnet_res   <- fread("/val_resnet.csv")
test_resnet_res  <- fread("/tma2_resnet.csv")

# =========================================================
# 3. Define constants and helper functions
# =========================================================
LEVELS <- paste0("CN", 1:10)

normalize_cn <- function(x) {
  x <- as.character(x)
  x <- ifelse(grepl("^CN", x), x, paste0("CN", x))
  sub("^CN0+", "CN", x)
}

make_confusion_matrix <- function(df) {
  true_label <- normalize_cn(df$label)
  pred_label <- normalize_cn(df$predicted_class)

  cm_df <- as.data.frame.array(table(true_label, pred_label))
  colnames(cm_df) <- paste0("CN", sub("^CN0*", "", colnames(cm_df)))
  rownames(cm_df) <- paste0("CN", sub("^CN0*", "", rownames(cm_df)))

  if ("CN0" %in% rownames(cm_df) && "CN0" %in% colnames(cm_df)) {
    cm_df <- cm_df[
      setdiff(rownames(cm_df), "CN0"),
      setdiff(colnames(cm_df), "CN0"),
      drop = FALSE
    ]
  }

  keep_rows <- intersect(rownames(cm_df), LEVELS)
  keep_cols <- intersect(colnames(cm_df), LEVELS)
  cm_df <- cm_df[keep_rows, keep_cols, drop = FALSE]

  cm <- matrix(0, nrow = 10, ncol = 10, dimnames = list(LEVELS, LEVELS))
  cm[rownames(cm_df), colnames(cm_df)] <- as.matrix(cm_df)
  cm
}

calc_accuracy_f1_kappa <- function(cm) {
  cm <- as.matrix(cm)
  total <- sum(cm)

  accuracy <- if (total == 0) NA_real_ else sum(diag(cm)) / total

  row_sums <- rowSums(cm)
  col_sums <- colSums(cm)
  valid_idx <- which(row_sums + col_sums > 0)

  calc_f1 <- function(i) {
    tp <- cm[i, i]
    fp <- sum(cm[, i]) - tp
    fn <- sum(cm[i, ]) - tp

    precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
    recall    <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)

    if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
      NA_real_
    } else {
      2 * precision * recall / (precision + recall)
    }
  }

  f1_values <- vapply(valid_idx, calc_f1, numeric(1))
  macro_f1 <- if (length(f1_values) == 0) NA_real_ else mean(f1_values, na.rm = TRUE)

  expected <- if (total == 0) NA_real_ else sum(row_sums * col_sums) / (total^2)
  kappa <- if (is.na(expected) || expected == 1) NA_real_ else (accuracy - expected) / (1 - expected)

  list(
    Accuracy = accuracy,
    MacroF1 = macro_f1,
    Kappa = kappa
  )
}

calc_multiclass_mcc <- function(cm) {
  cm <- as.matrix(cm)
  n_total <- sum(cm)
  if (n_total == 0) return(NA_real_)

  r <- rowSums(cm)
  c <- colSums(cm)

  numerator <- (sum(diag(cm)) * n_total) - sum(r * c)
  denominator <- sqrt((n_total^2 - sum(r^2)) * (n_total^2 - sum(c^2)))

  if (denominator == 0) return(NA_real_)
  numerator / denominator
}

detect_probability_matrix <- function(df) {
  nms <- names(df)

  get_prob <- function(fmt) {
    cols <- sprintf(fmt, LEVELS)
    if (all(cols %in% nms)) {
      as.matrix(df[, ..cols])
    } else {
      NULL
    }
  }

  prob_mat <- get_prob("prob_%s")
  if (!is.null(prob_mat)) return(prob_mat)

  prob_mat <- get_prob("p_%s")
  if (!is.null(prob_mat)) return(prob_mat)

  get_prob("%s_prob")
}

calc_macro_auroc <- function(df) {
  true_label <- normalize_cn(df$label)
  pred_label <- normalize_cn(df$predicted_class)

  keep <- (true_label %in% LEVELS) & (pred_label %in% LEVELS)
  true_label <- true_label[keep]
  pred_label <- pred_label[keep]

  if (length(true_label) == 0L) {
    return(NA_real_)
  }

  prob_mat <- detect_probability_matrix(df[keep, ])
  aucs <- c()

  if (!is.null(prob_mat)) {
    for (i in seq_along(LEVELS)) {
      y_true <- as.integer(true_label == LEVELS[i])
      if (length(unique(y_true)) < 2) next

      score <- prob_mat[, i]
      roc_obj <- try(pROC::roc(response = y_true, predictor = score, quiet = TRUE), silent = TRUE)

      if (!inherits(roc_obj, "try-error")) {
        aucs <- c(aucs, as.numeric(roc_obj$auc))
      }
    }
  } else {
    for (lv in LEVELS) {
      y_true <- as.integer(true_label == lv)
      if (length(unique(y_true)) < 2) next

      score <- as.integer(pred_label == lv)
      roc_obj <- try(pROC::roc(response = y_true, predictor = score, quiet = TRUE), silent = TRUE)

      if (!inherits(roc_obj, "try-error")) {
        aucs <- c(aucs, as.numeric(roc_obj$auc))
      }
    }
  }

  if (length(aucs) == 0) NA_real_ else mean(aucs, na.rm = TRUE)
}

evaluate_model <- function(df, model_name, set_name) {
  cm <- make_confusion_matrix(df)
  basic_metrics <- calc_accuracy_f1_kappa(cm)
  mcc <- calc_multiclass_mcc(cm)
  auroc <- calc_macro_auroc(df)

  data.frame(
    Set = set_name,
    Model = model_name,
    Accuracy = basic_metrics$Accuracy,
    MacroF1 = basic_metrics$MacroF1,
    Kappa = basic_metrics$Kappa,
    MCC = mcc,
    AUROC = auroc,
    row.names = NULL
  )
}

# =========================================================
# 4. Evaluate all models on validation and testing sets
# =========================================================
results <- rbind(
  evaluate_model(val_canvas_res,   "CANVAS",  "Validation"),
  evaluate_model(val_uni_res,      "uni",     "Validation"),
  evaluate_model(val_virchow_res,  "virchow", "Validation"),
  evaluate_model(val_resnet_res,   "resnet",  "Validation"),
  evaluate_model(test_canvas_res,  "CANVAS",  "Testing"),
  evaluate_model(test_uni_res,     "uni",     "Testing"),
  evaluate_model(test_virchow_res, "virchow", "Testing"),
  evaluate_model(test_resnet_res,  "resnet",  "Testing")
)

results$Model <- factor(results$Model, levels = c("CANVAS", "uni", "virchow", "resnet"))
results$Set   <- factor(results$Set, levels = c("Validation", "Testing"))

print(results)

# =========================================================
# 5. Export summary table
# =========================================================
fwrite(results, file = file.path(out_dir, "val.test_five.metrics_compare.csv"))

# =========================================================
# 6. Plot helper
# =========================================================
plot_metric <- function(df, metric, y_label, out_file, ylim = NULL) {
  p <- ggplot(df, aes(x = Set, y = .data[[metric]], fill = Model)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.9) +
    labs(x = NULL, y = y_label) +
    scale_fill_npg() +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )

  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }

  ggsave(filename = file.path(out_dir, out_file), plot = p, width = 5.6, height = 3.4, dpi = 300)
  p
}

# =========================================================
# 7. Generate metric plots
# =========================================================
accuracy_plot <- plot_metric(results, "Accuracy", "Accuracy", "metric_accuracy.pdf", ylim = c(0.7, 0.9))
macrof1_plot  <- plot_metric(results, "MacroF1", "Macro-F1", "metric_macroF1.pdf", ylim = c(0.7, 0.85))
kappa_plot    <- plot_metric(results, "Kappa", "Cohen's kappa", "metric_kappa.pdf", ylim = c(0.7, 0.85))
mcc_plot      <- plot_metric(results, "MCC", "MCC (multiclass)", "metric_mcc.pdf", ylim = c(0.7, 0.9))
auroc_plot    <- plot_metric(results, "AUROC", "AUROC (macro, OvA)", "metric_auroc.pdf", ylim = c(0.8, 0.95))

# =========================================================
# 8. Combine all plots
# =========================================================
combined_plot <- accuracy_plot + macrof1_plot + kappa_plot + mcc_plot + auroc_plot +
  plot_layout(ncol = 2)

ggsave(
  filename = file.path(out_dir, "all_metrics_combined.pdf"),
  plot = combined_plot,
  width = 9.25,
  height = 6.75,
  dpi = 300
)