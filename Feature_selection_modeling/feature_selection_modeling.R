library(Boruta)
library(caret)
library(corrr)
library(doMC)
library(doParallel)
library(dplyr)
library(foreach)
library(ggplot2)
library(glmnet)
library(gridExtra)
library(igraph)
library(parallel)
library(progressr)
library(randomForestSRC)
library(readr)
library(scales)
library(survcomp)
library(survival)
library(survminer)
library(tibble)
library(tidyverse)
library(timeROC)

save.path <- "/feat_selection_model_build/"
save.dir  <- paste0(save.path, "results/")
plot.dir  <- paste0(save.path, "plot/")

data <- read.csv(file = paste0(save.dir, "Lung_IO_feature_rawdata.csv"), header = TRUE, row.names = 1, check.names = FALSE)

feature.id <- colnames(data)[-c(1:3)]
feature.anno <- data.frame(
  feature.id = feature.id,
  anno = rep(NA, length(feature.id))
)

feature.anno$anno[1:10] <- "Composition"
feature.anno$anno[11:16] <- "Diversity"
feature.anno$anno[17:106] <- "Spatial_metric"
feature.anno$anno[107:206] <- "Interaction"
feature.anno$anno[207:261] <- "Distance"
feature.anno$anno[262] <- "Transition"

data <- data %>% filter(!is.na(pfs), !is.na(status))
data <- data[, colMeans(is.na(data)) < 0.3]
for (col in setdiff(colnames(data), c("pfs", "status"))) {
  data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)
}

x_raw <- data[, !(colnames(data) %in% c("pfs", "status"))]
x_raw <- x_raw[, sapply(x_raw, is.numeric)]
x <- scale(as.matrix(x_raw)) 
y <- Surv(data$pfs, data$status)

cor.mat <- cor(x, method = "spearman")
adj <- abs(cor.mat) > 0.95
g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
clust <- cluster_louvain(g)
rep.features <- sapply(unique(membership(clust)), function(cl) sample(names(membership(clust))[membership(clust) == cl], 1))
x_nocorr <- x[, rep.features]

###
##
#
feature_stats <- list()
for (feat in colnames(x_nocorr)) {
  marker <- x_nocorr[, feat]
  
  df_tmp <- data.frame(pfs = data$pfs, status = data$status, marker = marker)
  fit <- coxph(Surv(pfs, status) ~ marker, data = df_tmp)
  sfit <- summary(fit)
  
  hr <- sfit$coef[1, 2]
  hr_lower <- sfit$conf.int[,"lower .95"]
  hr_upper <- sfit$conf.int[,"upper .95"]
  pval <- sfit$coef[1, 5]
  
  lp <- predict(fit, type = "lp")
  
  cidx <- concordance.index(lp, df_tmp$pfs, df_tmp$status, outx = TRUE)
  
  auc_obj <- timeROC(T = df_tmp$pfs,
                     delta = df_tmp$status,
                     marker = lp,
                     cause = 1,
                     times = c(6, 12),
                     iid = TRUE)
  
  auc180 <- auc_obj$AUC[1]
  auc180.se <- auc_obj$inference$vect_sd_1[1]
  auc180.ci <- c(auc180 - 1.96 * auc180.se, auc180 + 1.96 * auc180.se)
  
  auc365 <- auc_obj$AUC[2]
  auc365.se <- auc_obj$inference$vect_sd_1[2]
  auc365.ci <- c(auc365 - 1.96 * auc365.se, auc365 + 1.96 * auc365.se)
  
  feature_stats[[feat]] <- data.frame(
    Feature = feat,
    HR = hr,
    HR.lower = hr_lower,
    HR.upper = hr_upper,
    P = pval,
    C_index = cidx$c.index,
    C_lower = cidx$lower,
    C_upper = cidx$upper,
    AUC_180 = auc180,
    AUC180_lower = auc180.ci[1],
    AUC180_upper = auc180.ci[2],
    AUC_365 = auc365,
    AUC365_lower = auc365.ci[1],
    AUC365_upper = auc365.ci[2]
  )
}

feature_stats.df <- bind_rows(feature_stats)
write.csv(feature_stats.df, paste0(save.dir, "all_feature_unicox_cindex_auc_summary.csv"), row.names = FALSE)

###
##
#
registerDoMC(cores = 20)

set.seed(123)
lambda_grid <- 10^seq(-4, -1.3, length = 50)
n_boot <- 100  
n_features <- ncol(x_nocorr)

lasso_freq <- rep(0, n_features)
names(lasso_freq) <- colnames(x_nocorr)

beta_matrix <- matrix(0, nrow = n_features, ncol = n_boot,
                      dimnames = list(colnames(x_nocorr), paste0("Iter", 1:n_boot)))

cvm_list <- list()  

pb <- txtProgressBar(min = 0, max = n_boot, style = 3)

# === Bootstrap LASSO-Cox ===
for (i in 1:n_boot) {
  idx <- sample(1:nrow(x_nocorr), replace = TRUE)
  
  fit <- cv.glmnet(
    x_nocorr[idx, ], y[idx],
    family = "cox", alpha = 1,
    lambda = lambda_grid, nfolds = 5,
    parallel = TRUE, standardize = FALSE,
    type.measure = "C"  
  )
  
  cvm_list[[i]] <- data.frame(
    Iteration = i,
    log_lambda = log(fit$lambda),
    C_index = fit$cvm
  )
  
  sel <- rownames(coef(fit, s = fit$lambda.1se))[coef(fit, s = fit$lambda.1se)[, 1] != 0]
  lasso_freq[sel] <- lasso_freq[sel] + 1
  
  coefs <- as.matrix(coef(fit, s = fit$lambda.1se))
  beta_matrix[names(coefs[,1]), i] <- coefs[,1]
  
  setTxtProgressBar(pb, i)
}
close(pb)

saveRDS(lasso_freq, file = paste0(save.dir, "lasso_freq.rds"))
saveRDS(beta_matrix, file = paste0(save.dir, "lasso_beta_matrix.rds"))


lasso_freq <- readRDS(paste0(save.dir, "lasso_freq.rds"))
beta_matrix <- readRDS(paste0(save.dir, "lasso_beta_matrix.rds"))


lasso_freq.df <- data.frame(
  term = names(lasso_freq),
  lasso_freq = lasso_freq,
  beta_mean = rowMeans(beta_matrix),
  beta_sd = apply(beta_matrix, 1, sd)
)

lasso_freq.df.od <- lasso_freq.df %>%
  arrange(desc(lasso_freq))
write.csv(lasso_freq.df.od, paste0(save.dir, "lasso_freq.df.od.csv"), row.names = FALSE)

cvm_df <- do.call(rbind, cvm_list)
write.csv(cvm_df, paste0(save.dir, "lasso_bootstrap_lambda_vs_cindex.csv"), row.names = FALSE)

###
##
#
set.seed(123)
rf_importance_list <- list()
rf_importance_df <- data.frame()

n_repeats <- 1000
pb_rf <- txtProgressBar(min = 0, max = n_repeats, style = 3)

for (i in 1:n_repeats) {
  idx <- sample(1:nrow(data), replace = TRUE)
  data_boot <- data[idx, c("pfs", "status", lasso_top)]
  data_boot$status <- as.factor(data_boot$status)
  
  rf_fit <- rfsrc(Surv(pfs, status) ~ ., data = data_boot, ntree = 500, importance = "permute")
  imp_vals <- rf_fit$importance[, "event.1"]
  imp_vals <- imp_vals[!is.na(imp_vals)]
  
  df_i <- data.frame(
    Feature = names(imp_vals),
    Importance = as.numeric(imp_vals),
    Iteration = i
  )
  rf_importance_df <- rbind(rf_importance_df, df_i)
  
  setTxtProgressBar(pb_rf, i)
}
close(pb_rf)

rf_importance_df#
write_excel_csv(rf_importance_df, file=paste0(save.dir, "times_1000_rf_importance_df.csv"))

save.dir  <- paste0(save.path, "results/")
plot.dir  <- paste0(save.path, "plot/")
rf_importance_df <- read.csv(file=paste0(save.dir, "times_1000_rf_importance_df.csv"), header = T, as.is = T)



format_p_scientific <- function(p) {
  if (is.na(p)) return(NA)
  if (p < 2e-16) {
    return("<2 × 10⁻¹⁶")
  } else {
    exp_str <- formatC(p, format = "e", digits = 1)  # 
    parts <- strsplit(exp_str, "e")[[1]]
    base <- parts[1]
    exp <- as.integer(parts[2])
    
    superscript <- c(
      `0` = "\u2070", `1` = "\u00B9", `2` = "\u00B2", `3` = "\u00B3",
      `4` = "\u2074", `5` = "\u2075", `6` = "\u2076", `7` = "\u2077",
      `8` = "\u2078", `9` = "\u2079", `-` = "\u207B"
    )
    exp_str_sup <- paste0(superscript[strsplit(as.character(exp), "")[[1]]], collapse = "")
    return(paste0(base, " × 10", exp_str_sup))
  }
}

###
##
#
format_p_scientific <- function(p) {
  if (is.na(p)) return(NA)
  if (p < 2e-16) return("< 2 × 10^-16")
  
  exp_str <- formatC(p, format = "e", digits = 1)
  parts <- strsplit(exp_str, "e")[[1]]
  base <- parts[1]
  exponent <- as.integer(parts[2])
  return(paste0(base, " × 10^", exponent))
}


rf_importance_summary <- rf_importance_df %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(meanImportance = mean(Importance, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(meanImportance))

top_features <- rf_importance_summary$Feature
rf_top_df <- rf_importance_df %>% filter(Feature %in% top_features)


imp_thresh <- quantile(rf_top_df$Importance, 0.25)  # Q1 

rf_freq_df <- rf_top_df %>%
  dplyr::filter(Importance > imp_thresh) %>%
  dplyr::count(Feature, name = "Freq") %>%
  tidyr::complete(Feature = top_features, fill = list(Freq = 0)) %>%
  { setNames(.$Freq, .$Feature) }


p_null <- mean(rf_freq_df) / n_repeats
rf_pval <- sapply(rf_freq_df, function(k) 1 - pbinom(k - 1, n_repeats, p_null))
rf_padj <- p.adjust(rf_pval, method = "bonferroni")


rf_top_summary <- rf_importance_summary %>%
  mutate(
    Freq = rf_freq_df[Feature],
    P_binomial = rf_pval[Feature],
    P_adj = rf_padj[Feature],
    P_label = sapply(P_binomial, format_p_scientific)
  )

rf_top_summary.df <- as.data.frame(rf_top_summary)
rf_top_summary.df.od <- rf_top_summary.df[order(rf_top_summary.df$P_binomial),]
rf_top_summary.df.od

write_excel_csv(rf_top_summary, file=paste0(save.dir, "times_1000_rf_top_summary.csv"))
rf_top_summary <- read.csv(file=paste0(save.dir, "times_1000_rf_top_summary.csv"), header = T, as.is = T)

###
##
#
library(dplyr)
library(survival)
library(scales)
library(survminer)

for (col in setdiff(colnames(data), c("pfs", "status"))) {
  data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)
}

selected_features <- c("cci_CN08_CN03", "cci_CN02_CN08", "cci_CN02_CN09", "F_mean_CN09", "Pair_corr_g_mean_CN02", "cci_CN04_CN09", "cci_CN01_CN03", "J_mean_CN03", "cci_CN09_CN03", "dis_CN04_CN08", "cci_CN05_CN03", "div_Richness", "frequency_CN06")

features <- intersect(selected_features, colnames(data))
final_formula <- as.formula(paste("Surv(pfs, status) ~", paste(features, collapse = " + ")))
final_model <- coxph(final_formula, data = data)

sum_model <- summary(final_model)
sum_model_comb <- cbind.data.frame(sum_model$coefficients, sum_model$conf.int[,-c(1:2)])
write.csv(sum_model_comb, file = paste0(save.dir, "final_cox_model_summary.csv"))

data$CoxScore_raw <- predict(final_model, type = "lp")
data$CoxScore <- rescale(data$CoxScore_raw, to = c(0, 1))  # 

data$RiskGroup <- ifelse(data$CoxScore > median(data$CoxScore), "1", "0")
data$RiskGroup <- factor(data$RiskGroup, levels = c("0", "1"))

fit <- survfit(Surv(pfs, status) ~ RiskGroup, data = data)

unicox <- coxph(Surv(pfs, status) ~ RiskGroup, data = data)
unicox_summary <- summary(unicox)

hr <- round(unicox_summary$coef[1, 2], 2)
lower <- round(unicox_summary$conf.int[,"lower .95"], 2)
upper <- round(unicox_summary$conf.int[,"upper .95"], 2)
pval <- signif(unicox_summary$coef[1, 5], 3)
pval_text <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", pval))
hr_text <- paste0("HR = ", hr, " (", lower, "–", upper, ")")

p_km <- ggsurvplot(
  fit, data = data,
  risk.table = TRUE,
  pval = TRUE,
  pval.coord = c(5, 0.05),
  palette = c("steelblue", "tomato"),
  title = "KM by Final Cox Risk Score Group",
  xlab = "Time (months)",
  ylab = "Survival probability",
  risk.table.height = 0.25,
  legend.labs = c("Low", "High"),
  legend.title = "Risk Group"
)

p_km$plot <- p_km$plot +
  annotate("text", x = 5, y = 0.15, label = hr_text, size = 3)

pdf(file = paste0(plot.dir, "KM_CoxRiskScore_Group.pdf"), width = 5.45, height = 6.89)
print(p_km)
dev.off()
