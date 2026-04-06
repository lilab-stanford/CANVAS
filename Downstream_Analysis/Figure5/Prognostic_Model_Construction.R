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