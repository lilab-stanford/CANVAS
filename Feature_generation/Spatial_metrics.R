library(data.table)
library(spatstat)
library(tidyr)
library(purrr)
library(dplyr)

#Spatial metrics: Aggregation and dispersion
save.path <- "rawdata/"
LUNG.IO_patch_level_rawdata <- fread(paste0(save.path, "LUNG.IO.rawdata.filter.final.csv"), header = T, check.names = F)
freq.rawdata <- LUNG.IO_patch_level_rawdata

freq.rawdata <- freq.rawdata %>%
  mutate(X_centroid = X_centroid / 1000,
         Y_centroid = Y_centroid / 1000)

calc_spatial_indices_safe <- function(data_sub) {
  coords <- data_sub[, c("X_centroid", "Y_centroid")]
  
  if (nrow(coords) < 3) {
    return(data.frame(
      Ripley_K_mean = NA, Ripley_L_mean = NA, Pair_corr_g_mean = NA, G_mean = NA,
      F_mean = NA, J_mean = NA, Clark_Evans = NA, Quadrat_chisq = NA, Kernel_density_mean = NA
    ))
  }
  
  win <- owin(xrange = range(coords$X_centroid), yrange = range(coords$Y_centroid))
  ppp_obj <- ppp(coords$X_centroid, coords$Y_centroid, window = win)
  
  r_max <- tryCatch({
    val <- max(bw.ppl(ppp_obj))
    if (!is.finite(val) || val <= 0) NA else val
  }, error = function(e) NA)
  
  if (is.na(r_max) || r_max <= 0) {
    r_max <- max(diff(range(coords$X_centroid)), diff(range(coords$Y_centroid))) / 10
  }
  
  if (is.na(r_max) || r_max <= 0) {
    r_max <- 1
  }
  
  min_spacing <- min(0.02, r_max / 500)
  r_vals <- seq(0, r_max, by = min_spacing)
  if (length(r_vals) < 2) {
    r_vals <- seq(0, 1, by = 0.01)
  }
  
  K_mean <- tryCatch({
    mean(Kest(ppp_obj, nlarge = 80000)$iso, na.rm = TRUE)
  }, error = function(e) NA)
  
  L_mean <- tryCatch({
    mean(Lest(ppp_obj, nlarge = 80000)$iso, na.rm = TRUE)
  }, error = function(e) NA)
  
  g_mean <- tryCatch({
    mean(pcf(ppp_obj, rmax = 20, divisor = "d")$iso, na.rm = TRUE)
  }, error = function(e) NA)
  
  G_mean <- tryCatch({
    mean(Gest(ppp_obj, r = r_vals)$rs, na.rm = TRUE)
  }, error = function(e) NA)
  
  F_mean <- tryCatch({
    mean(Fest(ppp_obj, r = r_vals)$rs, na.rm = TRUE)
  }, error = function(e) NA)
  
  J_mean <- tryCatch({
    mean(Jest(ppp_obj, r = r_vals)$rs, na.rm = TRUE)
  }, error = function(e) NA)
  
  clark_evans_index <- tryCatch({
    clarkevans(ppp_obj)[["Donnelly"]]
  }, error = function(e) NA)
  
  quadrat_chisq <- tryCatch({
    quadrat.test(ppp_obj, nx = 3, ny = 3)$statistic[[1]]
  }, error = function(e) NA)
  
  mean_density <- tryCatch({
    mean(density(ppp_obj)$v, na.rm = TRUE)
  }, error = function(e) NA)
  
  data.frame(
    Ripley_K_mean = K_mean,
    Ripley_L_mean = L_mean,
    Pair_corr_g_mean = g_mean,
    G_mean = G_mean,
    F_mean = F_mean,
    J_mean = J_mean,
    Clark_Evans = clark_evans_index,
    Quadrat_chisq = quadrat_chisq,
    Kernel_density_mean = mean_density
  )
}

all_results <- freq.rawdata %>%
  group_by(imageid, CN_label) %>%
  group_modify(~calc_spatial_indices_safe(.x)) %>%
  ungroup()

all_results

write.csv(all_results, file = paste0(save_path, "CN_spat_feature_results.csv"))


#chmod +x RCode.R
#nohup Rscript RCode.R > RCode.log 2>&1 &
