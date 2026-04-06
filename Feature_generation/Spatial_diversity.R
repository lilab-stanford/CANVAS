#diversity
library(dplyr)
library(vegan)
library(data.table)
library(tidyr)

save.path <- "rawdata/"
LUNG.IO_patch_level_rawdata <- fread(paste0(save.path, "LUNG.IO.rawdata.filter.final.csv"), header = T, check.names = F)
freq.rawdata <- LUNG.IO_patch_level_rawdata


diversity_indices <- freq.rawdata %>%
  group_by(imageid, CN_label) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = CN_label, values_from = count, values_fill = 0)

cn_cols <- colnames(diversity_indices)[!(colnames(diversity_indices) %in% "imageid")]

diversity_indices$Richness     <- vegan::specnumber(diversity_indices[, cn_cols])
diversity_indices$Shannon      <- vegan::diversity(diversity_indices[, cn_cols], index = "shannon")
diversity_indices$Simpson      <- vegan::diversity(diversity_indices[, cn_cols], index = "simpson")
diversity_indices$Inv_Simpson  <- vegan::diversity(diversity_indices[, cn_cols], index = "invsimpson")


diversity_indices$Fisher_alpha <- suppressWarnings(
  vegan::fisher.alpha(diversity_indices[, cn_cols])
)

diversity_indices$Pielou <- with(diversity_indices, Shannon / log(Richness + 1e-10))

final_diversity <- diversity_indices %>%
  select(imageid, Richness, Shannon, Simpson, Inv_Simpson, Pielou, Fisher_alpha)

print(summary(final_diversity))

dir.create(output.path, showWarnings = F)
write.csv(final_diversity, file = paste0(output.path, "CN_final_diversity_results.csv"), row.names = FALSE)


