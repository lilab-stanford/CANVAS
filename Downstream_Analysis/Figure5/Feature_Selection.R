rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(foreach)
  library(doMC)
  library(glmnet)
  library(igraph)
  library(randomForestSRC)
  library(readr)
  library(survival)
  library(survcomp)
  library(timeROC)
})

# =========================================================
# 1. Define paths
# =========================================================
save.path <- "/feat_selection_model_build/"
save.dir  <- file.path(save.path, "results")
plot.dir  <- file.path(save.path, "plot")

dir.create(save.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. Load data
# =========================================================
data <- read.csv(
  file = file.path(save.dir, "Lung_IO.csv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# =========================================================
# 3. Annotate feature groups
# =========================================================
feature.id <- colnames(data)[-c(1:3)]
feature.anno <- data.frame(
  feature.id = feature.id,
  anno = rep(NA_character_, length(feature.id)),
  stringsAsFactors = FALSE
)

feature.anno$anno[1:10]    <- "Composition"
feature.anno$anno[11:16]   <- "Diversity"
feature.anno$anno[17:106]  <- "Spatial_metric"
feature.anno$anno[107:206] <- "Interaction"
feature.anno$anno[207:261] <- "Distance"
feature.anno$anno[262]     <- "Transition"

write.csv(
  feature.anno,
  file = file.path(save.dir, "feature_annotation.csv"),
  row.names = FALSE
)

# =========================================================
# 4. Basic preprocessing
# =========================================================
data <- data %>% filter(!is.na(pfs), !is.na(status))
data <- data[, colMeans(is.na(data)) < 0.3, drop = FALSE]

for (col in setdiff(colnames(data), c("pfs", "status"))) {
  data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)
}

x_raw <- data[, !(colnames(data) %in% c("pfs", "status")), drop = FALSE]
x_raw <- x_raw[, sapply(x_raw, is.numeric), drop = FALSE]
x <- scale(as.matrix(x_raw))
y <- Surv(data$pfs, data$status)

# =========================================================
# 5. Remove highly correlated features
# =========================================================
set.seed(123)

cor.mat <- cor(x, method = "spearman")
adj <- abs(cor.mat) > 0.95
diag(adj) <- FALSE

g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
clust <- cluster_louvain(g)

rep.features <- sapply(
  unique(membership(clust)),
  function(cl) sample(names(membership(clust))[membership(clust) == cl], 1)
)

x_nocorr <- x[, rep.features, drop = FALSE]

write.csv(
  data.frame(feature = colnames(x_nocorr)),
  file = file.path(save.dir, "representative_features_after_correlation_filter.csv"),
  row.names = FALSE
)

# =========================================================
# 6. Univariate Cox, C-index, and time-dependent AUC
# =========================================================
feature_stats <- vector("list", length = ncol(x_nocorr))
names(feature_stats) <- colnames(x_nocorr)

for (feat in colnames(x_nocorr)) {
  marker <- x_nocorr[, feat]
  df_tmp <- data.frame(pfs = data$pfs, status = data$status, marker = marker)

  fit <- coxph(Surv(pfs, status) ~ marker, data = df_tmp)
  sfit <- summary(fit)

  hr <- sfit$coef[1, "exp(coef)"]
  hr_lower <- sfit$conf.int[1, "lower .95"]
  hr_upper <- sfit$conf.int[1, "upper .95"]
  pval <- sfit$coef[1, "Pr(>|z|)"]

  lp <- predict(fit, type = "lp")

  cidx <- concordance.index(
    x = lp,
    surv.time = df_tmp$pfs,
    surv.event = df_tmp$status,
    outx = TRUE
  )

  auc_obj <- timeROC(
    T = df_tmp$pfs,
    delta = df_tmp$status,
    marker = lp,
    cause = 1,
    times = c(6, 12),
    iid = TRUE
  )

  auc_6  <- auc_obj$AUC[1]
  auc_12 <- auc_obj$AUC[2]

  auc_6_se  <- auc_obj$inference$vect_sd_1[1]
  auc_12_se <- auc_obj$inference$vect_sd_1[2]

  feature_stats[[feat]] <- data.frame(
    Feature = feat,
    HR = hr,
    HR.lower = hr_lower,
    HR.upper = hr_upper,
    P = pval,
    C_index = cidx$c.index,
    C_lower = cidx$lower,
    C_upper = cidx$upper,
    AUC_6 = auc_6,
    AUC_6_lower = auc_6 - 1.96 * auc_6_se,
    AUC_6_upper = auc_6 + 1.96 * auc_6_se,
    AUC_12 = auc_12,
    AUC_12_lower = auc_12 - 1.96 * auc_12_se,
    AUC_12_upper = auc_12 + 1.96 * auc_12_se
  )
}

feature_stats.df <- bind_rows(feature_stats)

write.csv(
  feature_stats.df,
  file.path(save.dir, "all_feature_unicox_cindex_auc_summary.csv"),
  row.names = FALSE
)

# =========================================================
# 7. Bootstrap LASSO-Cox
# =========================================================
registerDoMC(cores = 20)

set.seed(123)
lambda_grid <- 10^seq(-4, -1.3, length = 50)
n_boot <- 100
n_features <- ncol(x_nocorr)

lasso_freq <- setNames(rep(0, n_features), colnames(x_nocorr))
beta_matrix <- matrix(
  0,
  nrow = n_features,
  ncol = n_boot,
  dimnames = list(colnames(x_nocorr), paste0("Iter", seq_len(n_boot)))
)

cvm_list <- vector("list", n_boot)
pb <- txtProgressBar(min = 0, max = n_boot, style = 3)

for (i in seq_len(n_boot)) {
  idx <- sample(seq_len(nrow(x_nocorr)), replace = TRUE)

  fit <- cv.glmnet(
    x_nocorr[idx, , drop = FALSE],
    y[idx],
    family = "cox",
    alpha = 1,
    lambda = lambda_grid,
    nfolds = 5,
    parallel = TRUE,
    standardize = FALSE,
    type.measure = "C"
  )

  cvm_list[[i]] <- data.frame(
    Iteration = i,
    log_lambda = log(fit$lambda),
    C_index = fit$cvm
  )

  coef_mat <- as.matrix(coef(fit, s = fit$lambda.1se))
  selected <- rownames(coef_mat)[coef_mat[, 1] != 0]

  if (length(selected) > 0) {
    lasso_freq[selected] <- lasso_freq[selected] + 1
  }

  beta_matrix[rownames(coef_mat), i] <- coef_mat[, 1]
  setTxtProgressBar(pb, i)
}
close(pb)

saveRDS(lasso_freq, file = file.path(save.dir, "lasso_freq.rds"))
saveRDS(beta_matrix, file = file.path(save.dir, "lasso_beta_matrix.rds"))

lasso_freq.df <- data.frame(
  Feature = names(lasso_freq),
  LASSO_Frequency = as.numeric(lasso_freq),
  Beta_Mean = rowMeans(beta_matrix),
  Beta_SD = apply(beta_matrix, 1, sd),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(LASSO_Frequency))

write.csv(
  lasso_freq.df,
  file.path(save.dir, "lasso_freq_summary.csv"),
  row.names = FALSE
)

cvm_df <- bind_rows(cvm_list)
write.csv(
  cvm_df,
  file.path(save.dir, "lasso_bootstrap_lambda_vs_cindex.csv"),
  row.names = FALSE
)

# =========================================================
# 8. Select top LASSO features for RF survival importance
# =========================================================
lasso_top <- lasso_freq.df %>%
  filter(LASSO_Frequency > 0) %>%
  arrange(desc(LASSO_Frequency)) %>%
  pull(Feature)

write.csv(
  data.frame(Feature = lasso_top),
  file.path(save.dir, "lasso_selected_features.csv"),
  row.names = FALSE
)

# =========================================================
# 9. Repeated random survival forest importance
# =========================================================
set.seed(123)

rf_importance_df <- data.frame()
n_repeats <- 1000
pb_rf <- txtProgressBar(min = 0, max = n_repeats, style = 3)

for (i in seq_len(n_repeats)) {
  idx <- sample(seq_len(nrow(data)), replace = TRUE)

  data_boot <- data[idx, c("pfs", "status", lasso_top), drop = FALSE]

  rf_fit <- rfsrc(
    Surv(pfs, status) ~ .,
    data = data_boot,
    ntree = 500,
    importance = "permute"
  )

  imp_vals <- rf_fit$importance[, "event.1"]
  imp_vals <- imp_vals[!is.na(imp_vals)]

  df_i <- data.frame(
    Feature = names(imp_vals),
    Importance = as.numeric(imp_vals),
    Iteration = i,
    stringsAsFactors = FALSE
  )

  rf_importance_df <- rbind(rf_importance_df, df_i)
  setTxtProgressBar(pb_rf, i)
}
close(pb_rf)

write_excel_csv(
  rf_importance_df,
  file = file.path(save.dir, "rf_importance_1000_iterations.csv")
)

# =========================================================
# 10. Summarize RF importance and significance
# =========================================================
format_p_scientific <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 2e-16) return("< 2 × 10^-16")

  exp_str <- formatC(p, format = "e", digits = 1)
  parts <- strsplit(exp_str, "e")[[1]]
  base <- parts[1]
  exponent <- as.integer(parts[2])

  paste0(base, " × 10^", exponent)
}

rf_importance_summary <- rf_importance_df %>%
  group_by(Feature) %>%
  summarise(
    Mean_Importance = mean(Importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Importance))

top_features <- rf_importance_summary$Feature
rf_top_df <- rf_importance_df %>% filter(Feature %in% top_features)

imp_thresh <- quantile(rf_top_df$Importance, 0.25, na.rm = TRUE)

rf_freq_df <- rf_top_df %>%
  filter(Importance > imp_thresh) %>%
  count(Feature, name = "Freq") %>%
  complete(Feature = top_features, fill = list(Freq = 0))

p_null <- mean(rf_freq_df$Freq) / n_repeats
rf_pval <- sapply(rf_freq_df$Freq, function(k) 1 - pbinom(k - 1, n_repeats, p_null))
rf_padj <- p.adjust(rf_pval, method = "bonferroni")

rf_top_summary <- rf_importance_summary %>%
  left_join(rf_freq_df, by = "Feature") %>%
  mutate(
    P_binomial = rf_pval[match(Feature, rf_freq_df$Feature)],
    P_adj = rf_padj[match(Feature, rf_freq_df$Feature)],
    P_label = sapply(P_binomial, format_p_scientific)
  ) %>%
  arrange(P_binomial)

write_excel_csv(
  rf_top_summary,
  file = file.path(save.dir, "rf_top_summary_1000_iterations.csv")
)