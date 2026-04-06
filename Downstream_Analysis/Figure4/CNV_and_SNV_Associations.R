rm(list = ls())

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ordinal)
  library(emmeans)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(reshape2)
})

# =========================================================
# 1. Define paths
# =========================================================
base_dir <- "/results"
ccp_file <- file.path(base_dir, "/tcga.luad_sam.inf_Cluster_res.csv")

cnv_raw_file <- "/TCGA-LUAD_CNV.txt"
gene_anno_file <- "E:/lyc.data/gene.id.annotation/gencode_v22_genetype.csv"
cnv_out_dir <- file.path(base_dir, "tcga.luad.cnv/results")
cnv_plot_dir <- file.path(base_dir, "tcga.luad.cnv/plot")

snv_raw_dir <- file.path(base_dir, "tcga.luad.mutation/rawdata")
snv_raw_file <- file.path(snv_raw_dir, "LUAD_mc3_gene_level.csv")
snv_out_dir <- file.path(base_dir, "tcga.luad.mutation/results")
snv_plot_dir <- file.path(base_dir, "tcga.luad.mutation/plot")

pan_driver_file <- "E:/F.new/文献/pan-cancer/Comprehensive Characterization of Cancer Driver Genes and Mutations/1-s2.0-S009286741830237X-mmc1.csv"
gene_pathway_anno_file <- "E:/F.new/文献/pan-cancer/Comprehensive Characterization of Cancer Driver Genes and Mutations/gene.sig.anno.csv"

dir.create(cnv_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cnv_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(snv_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(snv_plot_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. Optional: download TCGA LUAD mutation data
# =========================================================
# query <- GDCquery(
#   project = "TCGA-LUAD",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# setwd(snv_raw_dir)
# GDCdownload(query)
# GDCprepare(query, save = TRUE, save.filename = "TCGA-LUAD_SNP.Rdata")

# =========================================================
# 3. Load cluster assignment
# =========================================================
tcga.luad_sam.inf_Cluster_res <- read.csv(
  file = ccp_file,
  header = TRUE,
  as.is = TRUE
)

# =========================================================
# 4. Load and prepare CNV matrix
# =========================================================
tcga.luad.cnv.res <- read.table(
  file = cnv_raw_file,
  sep = "\t",
  check.names = FALSE,
  header = TRUE
)

# =========================================================
# 5. Summarize CNV frequency
# =========================================================
data <- tcga.luad.cnv.res.1

freq.loss <- apply(data[, 4:ncol(data)], 1, function(x) sum(x == -1, na.rm = TRUE))
freq.gain <- apply(data[, 4:ncol(data)], 1, function(x) sum(x == 1, na.rm = TRUE))

freq.results <- data.frame(
  Gene_ID = data$`Gene ID`,
  Loss_Frequency = freq.loss,
  Gain_Frequency = freq.gain
)

freq.results$freq <- (freq.results$Loss_Frequency + freq.results$Gain_Frequency) / length(inter.sam.id)
write.csv(freq.results, file = file.path(cnv_out_dir, "tcga.luad_cnv_frequency_summary.csv"), row.names = FALSE)

# =========================================================
# 6. Prepare long-format CNV table
# =========================================================
sample_info <- tcga.luad_sam.inf_Cluster_res.1
colnames(sample_info)[1] <- "sample_id"

long_data <- pivot_longer(
  data,
  cols = starts_with("TCGA"),
  names_to = "sample_id",
  values_to = "cnv_value"
)

long_data <- left_join(long_data, sample_info, by = "sample_id")

# =========================================================
# 7. Cluster-associated CNV analysis using ordinal model
# =========================================================
final_results <- data.frame(
  Gene_ID = character(),
  CNV_Type = integer(),
  Cluster = character(),
  Odds_Ratio = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (gene in unique(long_data$`Gene ID`)) {
  gene_data <- long_data[long_data$`Gene ID` == gene, ]
  gene_data <- gene_data[gene_data$cnv_value %in% c(-1, 1), ]

  if (nrow(gene_data) == 0) next
  if (nlevels(factor(gene_data$cnv_value)) < 2) next

  gene_data$cnv_value <- factor(gene_data$cnv_value, levels = c(-1, 1))

  model <- tryCatch(
    clm(cnv_value ~ factor(consensusClass), data = gene_data, link = "logit"),
    error = function(e) NULL
  )
  if (is.null(model)) next

  emm_res <- tryCatch(
    emmeans(model, specs = ~ consensusClass, adjust = "none"),
    error = function(e) NULL
  )
  if (is.null(emm_res)) next

  summary_emm <- summary(emm_res, type = "response", infer = c(TRUE, TRUE))

  for (i in seq_len(nrow(summary_emm))) {
    final_results <- rbind(
      final_results,
      data.frame(
        Gene_ID = gene,
        CNV_Type = NA_integer_,
        Cluster = summary_emm$consensusClass[i],
        Odds_Ratio = exp(summary_emm$emmean[i]),
        Lower_CI = exp(summary_emm$asymp.LCL[i]),
        Upper_CI = exp(summary_emm$asymp.UCL[i]),
        pvalue = summary_emm$p.value[i]
      )
    )
  }
}

final_results <- cbind.data.frame(
  tcga.luad.cnv.res.1[
    match(final_results$Gene_ID, tcga.luad.cnv.res.1$`Gene ID`),
    1:3
  ],
  final_results
)

write.csv(
  final_results,
  file = file.path(cnv_out_dir, "tcga.luad_cluster.diff.res.csv"),
  row.names = FALSE
)

# =========================================================
# 8. Select genes with elevated CNV odds ratio
# =========================================================
cnv_or_data <- final_results[!is.na(final_results$pvalue), ]

selected_genes <- cnv_or_data %>%
  group_by(Gene.ID) %>%
  summarise(Max_OR = max(Odds_Ratio, na.rm = TRUE), .groups = "drop") %>%
  filter(Max_OR > 1)

result <- cnv_or_data %>%
  filter(Gene.ID %in% selected_genes$Gene.ID) %>%
  na.omit()

write.csv(result, file = file.path(cnv_out_dir, "result.csv"), row.names = FALSE)

# =========================================================
# 9. KEGG enrichment for selected CNV genes
# =========================================================
gene_list <- unique(result$Gene.ID)

geneID <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = TRUE
)

if (!is.null(geneID) && nrow(geneID) > 0) {
  kegg_result <- enrichKEGG(
    gene = geneID$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05
  )

  eKEGG <- setReadable(kegg_result, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  diff.genes.kegg.res <- eKEGG@result

  write.csv(
    diff.genes.kegg.res,
    file = file.path(cnv_out_dir, "tcga.luad.CNV.OR.diff.genes.kegg.res.csv"),
    row.names = FALSE
  )

  s.sig.path <- c(
    "Estrogen signaling pathway",
    "mTOR signaling pathway",
    "Rap1 signaling pathway",
    "Cellular senescence",
    "Apoptosis - multiple species",
    "Hippo signaling pathway",
    "Drug metabolism - cytochrome P450",
    "EGFR tyrosine kinase inhibitor resistance",
    "HIF-1 signaling pathway",
    "TNF signaling pathway",
    "Steroid biosynthesis",
    "Signaling pathways regulating pluripotency of stem cells"
  )

  luad.related.diff.genes.kegg.res <- diff.genes.kegg.res[
    match(s.sig.path, diff.genes.kegg.res$Description),
  ]

  write.csv(
    luad.related.diff.genes.kegg.res,
    file = file.path(cnv_out_dir, "luad.related.diff.genes.kegg.res.csv"),
    row.names = FALSE
  )
}

# =========================================================
# 10. Define pathway-focused CNV gene panel
# =========================================================
mTOR_signaling_pathway.id <- c("MTOR", "AKT3", "PIK3CA")
Rap1_signaling_pathway.id <- c("KRAS", "RAP1GAP")
Cellular_senescence.id <- c("CDKN1A")
Apoptosis__multiple_species.id <- c("BCL2", "CASP8")
Hippo_signaling_pathway.id <- c("YAP1", "WWC1")
Drug_metabolism__cytochrome_P450.id <- c("CYP3A4")
EGFR_tyrosine_kinase_inhibitor_resistance.id <- c("EGFR", "ERBB2", "BRAF")
HIF1_signaling_pathway.id <- c("HIF1A")
TNF_signaling_pathway.id <- c("TNFRSF1A", "TRAF5")
Steroid_biosynthesis.id <- c("CYP51A1")
Signaling_pathways_regulating_pluripotency_of_stem_cells.id <- c("NANOG", "SOX2")
Immuno_checkpoint.id <- c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "TIGIT")
Other.id <- c("TGFB2")

all.path.gene.id <- c(
  mTOR_signaling_pathway.id,
  Rap1_signaling_pathway.id,
  Cellular_senescence.id,
  Apoptosis__multiple_species.id,
  Hippo_signaling_pathway.id,
  Drug_metabolism__cytochrome_P450.id,
  EGFR_tyrosine_kinase_inhibitor_resistance.id,
  HIF1_signaling_pathway.id,
  TNF_signaling_pathway.id,
  Steroid_biosynthesis.id,
  Signaling_pathways_regulating_pluripotency_of_stem_cells.id,
  Immuno_checkpoint.id,
  Other.id
)

all.path.gene.id.OR.result <- result[result$Gene.ID %in% all.path.gene.id, ]
sorted_result <- all.path.gene.id.OR.result[
  order(match(all.path.gene.id.OR.result$Gene.ID, all.path.gene.id)),
]

write.csv(
  sorted_result,
  file = file.path(cnv_out_dir, "tcga.luad.CNV.OR.all.path.gene.id.comb.res.csv"),
  row.names = FALSE
)

# =========================================================
# 11. Plot CNV log2 odds ratio
# =========================================================

p_cnv <- ggplot(g.data, aes(x = Gene.ID, y = lg2.OR, fill = color.group)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_identity() +
  facet_wrap(~ Cluster, scales = "free_x", strip.position = "top", nrow = 1, ncol = 4) +
  labs(
    title = "Log2 Odds Ratio by Gene",
    x = "Gene",
    y = "Log2 Odds Ratio"
  ) +
  theme_classic()

pdf(file.path(cnv_plot_dir, "each.cluster_CNV.gene.sig_plot.pdf"), width = 5.09, height = 6.85)
print(p_cnv)
dev.off()

# =========================================================
# 12. Load and prepare SNV matrix
# =========================================================
tcga.luad.mutation.res <- read.csv(
  file = snv_raw_file,
  header = TRUE,
  sep = ",",
  check.names = FALSE
)

tcga.luad.mutation.res <- tcga.luad.mutation.res[!duplicated(tcga.luad.mutation.res$sample), ]
tcga.luad.mutation.res <- na.omit(tcga.luad.mutation.res)
rownames(tcga.luad.mutation.res) <- tcga.luad.mutation.res$sample

colnames(tcga.luad.mutation.res) <- substring(colnames(tcga.luad.mutation.res), 1, 12)
colnames(tcga.luad.mutation.res)[1] <- "Gene.Symbol"

inter.sam.id <- intersect(tcga.luad_sam.inf_Cluster_res$X, colnames(tcga.luad.mutation.res))

tcga.luad_sam.inf_Cluster_res.1 <- tcga.luad_sam.inf_Cluster_res[
  match(inter.sam.id, tcga.luad_sam.inf_Cluster_res$X),
]

tcga.luad.mutation.res <- cbind.data.frame(
  Gene.Symbol = tcga.luad.mutation.res$Gene.Symbol,
  tcga.luad.mutation.res[, inter.sam.id, drop = FALSE]
)

# =========================================================
# 13. Summarize mutation frequency by cluster
# =========================================================
cluster_ids <- list(
  Cluster_1 = tcga.luad_sam.inf_Cluster_res.1[tcga.luad_sam.inf_Cluster_res.1$consensusClass == "Cluster_1", "bcr_patient_barcode"],
  Cluster_2 = tcga.luad_sam.inf_Cluster_res.1[tcga.luad_sam.inf_Cluster_res.1$consensusClass == "Cluster_2", "bcr_patient_barcode"],
  Cluster_3 = tcga.luad_sam.inf_Cluster_res.1[tcga.luad_sam.inf_Cluster_res.1$consensusClass == "Cluster_3", "bcr_patient_barcode"],
  Cluster_4 = tcga.luad_sam.inf_Cluster_res.1[tcga.luad_sam.inf_Cluster_res.1$consensusClass == "Cluster_4", "bcr_patient_barcode"]
)

long_mutation_data <- pivot_longer(
  tcga.luad.mutation.res,
  cols = -Gene.Symbol,
  names_to = "sample_id",
  values_to = "Mutation_Status"
)

long_mutation_data$Mutation_Status <- as.numeric(long_mutation_data$Mutation_Status)

gene_frequencies <- lapply(names(cluster_ids), function(cluster) {
  sample_ids <- cluster_ids[[cluster]]

  long_mutation_data %>%
    filter(sample_id %in% sample_ids) %>%
    group_by(Gene.Symbol) %>%
    summarize(
      Mutation_Count = sum(Mutation_Status, na.rm = TRUE),
      Total_Samples = n(),
      Frequency = Mutation_Count / Total_Samples * 100,
      .groups = "drop"
    ) %>%
    mutate(Cluster = cluster)
})

cluster.4.mutation.freq.final_results <- bind_rows(gene_frequencies)

write.csv(
  cluster.4.mutation.freq.final_results,
  file = file.path(snv_out_dir, "cluster.4.mutation.freq.final_results.csv"),
  row.names = FALSE
)

# =========================================================
# 14. Select mutation genes of interest
# =========================================================
pan.cancer.mutation.gene.id <- read.csv(
  file = pan_driver_file,
  header = TRUE,
  as.is = TRUE
)

all.path.gene.id <- unique(pan.cancer.mutation.gene.id$Gene)

mutation.diff.s.gene.res <- as.data.frame(
  cluster.4.mutation.freq.final_results[
    cluster.4.mutation.freq.final_results$Gene.Symbol %in% all.path.gene.id,
  ]
)

mutation.diff.s.gene.res <- mutation.diff.s.gene.res[order(mutation.diff.s.gene.res$Gene.Symbol), ]

# =========================================================
# 15. Fisher test for cluster-associated mutation genes
# =========================================================
CMH.rawdata <- cbind.data.frame(
  consensusClass = tcga.luad_sam.inf_Cluster_res.1[, 2],
  t(tcga.luad.mutation.res[, -1])
)

gene_freq <- colSums(CMH.rawdata[, 2:ncol(CMH.rawdata)], na.rm = TRUE) / nrow(CMH.rawdata)
selected_genes <- names(gene_freq[gene_freq >= 0.01])

filtered_data <- CMH.rawdata[, c("consensusClass", selected_genes), drop = FALSE]

results <- data.frame(Gene = character(), pvalue = numeric(), stringsAsFactors = FALSE)

genes <- colnames(filtered_data)[2:ncol(filtered_data)]

for (gene in genes) {
  table_data <- table(filtered_data$consensusClass, filtered_data[[gene]])

  if (all(table_data == 0)) {
    results <- rbind(results, data.frame(Gene = gene, pvalue = 1))
    next
  }

  fisher_result <- tryCatch(
    fisher.test(table_data),
    error = function(e) list(p.value = 1)
  )

  results <- rbind(
    results,
    data.frame(Gene = gene, pvalue = fisher_result$p.value)
  )
}

results$FDR <- p.adjust(results$pvalue, method = "BH")

fisher.diff.path.gene.results <- results[match(all.path.gene.id, results$Gene), ]
fisher.diff.path.gene.results <- na.omit(fisher.diff.path.gene.results)
fisher.diff.path.gene.results.sig <- fisher.diff.path.gene.results[
  fisher.diff.path.gene.results$pvalue < 0.05,
]

write.csv(
  fisher.diff.path.gene.results.sig,
  file = file.path(snv_out_dir, "fisher.diff.path.gene.results.sig.csv"),
  row.names = FALSE
)

# =========================================================
# 16. Build mutation heatmap input
# =========================================================
cluster.4.mutation.freq.sig.gene.res <- cluster.4.mutation.freq.final_results[
  cluster.4.mutation.freq.final_results$Gene.Symbol %in% fisher.diff.path.gene.results.sig$Gene,
]

cluster.4.mutation.freq.sig.gene.res <- cluster.4.mutation.freq.sig.gene.res[
  order(cluster.4.mutation.freq.sig.gene.res$Gene.Symbol),
]

g.data <- cluster.4.mutation.freq.sig.gene.res
g.data.wide <- dcast(g.data, Gene.Symbol ~ Cluster, value.var = "Frequency")

gene.sig.anno <- read.csv(
  file = gene_pathway_anno_file,
  header = TRUE,
  as.is = TRUE
)

g.data.wide <- cbind.data.frame(
  Pathway = gene.sig.anno[match(g.data.wide$Gene.Symbol, gene.sig.anno$Gene), ][, 3],
  g.data.wide
)
g.data.wide$Pathway[is.na(g.data.wide$Pathway)] <- "Others pathways"

gene.anno <- g.data.wide
write.csv(gene.anno, file = file.path(snv_out_dir, "gene.anno.csv"), row.names = FALSE)

# =========================================================
# 17. Plot mutation frequency heatmap
# =========================================================
g.data.wide <- g.data.wide[, -1, drop = FALSE]
g.data.wide[, -1] <- t(apply(g.data.wide[, -1, drop = FALSE], 1, scale))

g.data.wide$Gene.Symbol <- factor(g.data.wide$Gene.Symbol, levels = rev(g.data.wide$Gene.Symbol))

g.data.long <- melt(
  g.data.wide,
  id.vars = "Gene.Symbol",
  variable.name = "Cluster",
  value.name = "Normalized_Frequency"
)

p_snv <- ggplot(g.data.long, aes(x = Cluster, y = Gene.Symbol, fill = Normalized_Frequency)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#4F81BD",
    high = "#C0504D",
    mid = "white",
    midpoint = 0,
    limit = c(
      min(g.data.long$Normalized_Frequency, na.rm = TRUE),
      max(g.data.long$Normalized_Frequency, na.rm = TRUE)
    ),
    space = "Lab",
    name = "Z-Score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 12, face = "bold")
  )

pdf(file.path(snv_plot_dir, "each.cluster_mut.gene.freq_plot.pdf"), width = 3.08, height = 6.85)
print(p_snv)
dev.off()