library(dplyr)
library(tidyr)
library(purrr)

###
##
#
save.path <- "rawdata/"
LUNG.IO_patch_level_rawdata <- fread(paste0(save.path, "LUNG.IO.rawdata.filter.final.csv"), header = T, check.names = F)
freq.rawdata <- LUNG.IO_patch_level_rawdata

cn_summary <- freq.rawdata %>%
  group_by(imageid, CN_label) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(imageid) %>%
  mutate(frequency = count / sum(count),
         total_count = sum(count)) %>%
  pivot_wider(names_from = CN_label,
              values_from = c(count, frequency),
              names_sep = "_",
              values_fill = 0)

cn_labels <- c("CN01", "CN04", "CN06", "CN02", "CN03", "CN05", "CN08", "CN09", "CN10", "CN07")

cn_combinations_2 <- combn(cn_labels, 1, simplify = FALSE)


all_combinations <- c(cn_combinations_2)

cn_combination_freq <- freq.rawdata %>%
  group_by(imageid) %>%
  summarise(total_count = n(), .groups = 'drop') %>%
  bind_cols(map_dfc(all_combinations, function(comb) {
    freq.rawdata %>%
      group_by(imageid) %>%
      summarise(freq = sum(CN_label %in% comb) / n(), .groups = 'drop') %>%
      select(freq) %>%
      setNames(paste0("freq_", paste(comb, collapse = "_")))
  }))
cn_combination_freq <- cn_combination_freq[,-which(colnames(cn_combination_freq)=="total_count")]
final_cn_summary <- left_join(cn_summary, cn_combination_freq, by = "imageid")

head(final_cn_summary)

write.csv(final_cn_summary, file=paste0(save_path, "LUNG.IO_count.freq.res.csv"))
