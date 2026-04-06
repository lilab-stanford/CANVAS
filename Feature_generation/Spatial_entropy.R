library(data.table)
library(dplyr)
library(FNN)
library(entropy)
library(spatstat.geom)
library(spatstat.core)
library(SpatEntropy)

save.path <- "rawdata/"
LUNG.IO_patch_level_rawdata <- fread(paste0(save.path, "LUNG.IO.rawdata.filter.final.csv"), header = T, check.names = F)
df.all <- LUNG.IO_patch_level_rawdata

entropy_results <- data.frame()

image.ids <- unique(df.all$imageid)

for (img in image.ids) {
  cat("Processing image:", img, "\n")
  
  df <- df.all %>% filter(imageid == img)
  
  result_row <- data.frame(
    imageid = img,
    ShannonZ = 0,
    Altieri_SMI = 0,
    Altieri_RES = 0,
    Leibovici = 0
  )
  
  if (nrow(df) < 10 || length(unique(df$CN_label)) < 2) {
    entropy_results <- rbind(entropy_results, result_row)
    next
  }
  
  df$CN_label <- as.factor(df$CN_label)
  xrange <- range(df$X_centroid)
  yrange <- range(df$Y_centroid)
  win <- owin(xrange = xrange, yrange = yrange)
  pp <- ppp(x = df$X_centroid, y = df$Y_centroid, window = win, marks = df$CN_label)
  
  mindist <- tryCatch(quantile(nndist(pp), probs = 0.1), error = function(e) NA)
  if (is.na(mindist) || mindist == 0) {
    entropy_results <- rbind(entropy_results, result_row)
    next
  }
  distbreak <- c(mindist, mindist * 2)
  
  alt <- tryCatch(altieri(pp, distbreak = distbreak, plotout = FALSE), error = function(e) NULL)
  leib <- tryCatch(leibovici(pp, mindist), error = function(e) NULL)
  
  if (!is.null(alt) && !is.null(leib) &&
      length(alt$ShannonZ$shannZ) == 1 &&
      length(alt$SMI) == 1 &&
      length(alt$RES) == 1 &&
      length(leib$entropy) == 1) {
    
    result_row$ShannonZ <- alt$ShannonZ$shannZ
    result_row$Altieri_SMI <- alt$SMI
    result_row$Altieri_RES <- alt$RES
    result_row$Leibovici <- leib$entropy
    cat("Image", img, "entropy computed.\n")
  } else {
    cat("Image", img, "entropy computation failed. Set to 0.\n")
  }
  
  entropy_results <- rbind(entropy_results, result_row)
}

write.csv(entropy_results, paste0(save.path, "spatial_entropy_summary_single_thread.csv"), row.names = FALSE)
cat("All images processed and saved.\n")


###
##
#
save.path <- "rawdata/"
LUNG.IO_patch_level_rawdata <- fread(paste0(save.path, "LUNG.IO.rawdata.filter.final.csv"), header = T, check.names = F)
df.all <- LUNG.IO_patch_level_rawdata

compute_spatial_transition_entropy <- function(df.sub, k = 6) {
  if (nrow(df.sub) < k + 1 || length(unique(df.sub$CN_label)) < 2) return(0)
  
  coords <- cbind(df.sub$X_centroid, df.sub$Y_centroid)
  labels <- as.character(df.sub$CN_label)
  n <- nrow(df.sub)
  
  nn <- get.knn(coords, k = k)$nn.index
  
  label_set <- sort(unique(labels))
  transitions <- matrix(0, nrow = length(label_set), ncol = length(label_set),
                        dimnames = list(label_set, label_set))
  
  for (i in 1:n) {
    from <- labels[i]
    to <- labels[nn[i, ]]
    for (t in to) {
      transitions[from, t] <- transitions[from, t] + 1
    }
  }
  
  P <- transitions / rowSums(transitions + 1e-10)
  entropy_matrix <- -P * log(P + 1e-10)
  ste <- sum(entropy_matrix, na.rm = TRUE)
  return(ste)
}

image.ids <- unique(df.all$imageid)
ste.results <- data.frame()

for (img in image.ids) {
  cat("Processing image:", img, "\n")
  df.img <- df.all %>% filter(imageid == img)
  ste <- compute_spatial_transition_entropy(df.img, k = 6)
  ste.results <- rbind(ste.results, data.frame(imageid = img, SpatialTransitionEntropy = ste))
}

write.csv(ste.results, paste0(save.path, "spatial_transition_entropy.csv"), row.names = FALSE)
cat("All STE values computed and saved.\n")

