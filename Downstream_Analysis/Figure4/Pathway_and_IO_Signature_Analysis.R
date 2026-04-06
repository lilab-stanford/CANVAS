rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(GSVA)
  library(GSEABase)
  library(clusterProfiler)
  library(limma)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(mclust)
  library(DescTools)
  library(ggalluvial)
  library(car)
})

# =========================================================
# 1. Define paths
# =========================================================
base_dir <- "/cluster_subtype/scimap/results/luad_subtype_omics"

raw_dir     <- file.path(base_dir, "rawdata")
result_dir  <- file.path(base_dir, "results")
plot_dir    <- file.path(result_dir, "plot")
gsea_dir    <- file.path(result_dir, "GESA_4.diff.res")
gsea_plot_dir <- file.path(gsea_dir, "plot")
cyto_dir    <- file.path(result_dir, "cytoscape")
sankey_dir  <- file.path(base_dir, "sankey.plot")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cyto_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sankey_dir, recursive = TRUE, showWarnings = FALSE)

# External inputs
io_sig_file <- "/IO_subtype_gene.signature_Cancer.cell.csv"
io_gmt_file <- "/IO.sig_gene.set_Cancell.gmt"
tpm_file    <- "F:/TCGA_GTEx/TCGA_file/TPM_all/TCGA_LUAD_TPM.txt"
gene_anno_file <- "E://gene.id.annotation/gencode_v22_genetype.csv"
cluster_file <- "/results/tcga.luad_sam.inf_Cluster_res.csv"
io_subtype_inf_file <- "/TCGA_IO_subtype_inf.csv"
io_quant_score_file <- "/TCGA_pan.cancer_Geno.Trans_IO.signature_comb.csv"

hallmark_gmt <- "E:/GSEA.signature/2025/h.all.v2024.1.Hs.symbols.gmt"
kegg_gmt     <- "E:/GSEA.signature/2025/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"
reactome_gmt <- "E:/GSEA.signature/2025/c2.cp.reactome.v2024.1.Hs.symbols.gmt"
biocarta_gmt <- "E:/GSEA.signature/2025/c2.cp.biocarta.v2024.1.Hs.symbols.gmt"

# =========================================================
# 2. Build IO gene-set GMT file
# =========================================================
io.gene.sig <- read.csv(io_sig_file, header = TRUE, as.is = TRUE)
tmp <- io.gene.sig[, c(2, 1)]
colnames(tmp) <- c("tissue", "gene.id")

io_gene_list <- tapply(tmp$gene.id, as.factor(tmp$tissue), function(x) unique(x))

write_gmt <- function(gene_set, gmt_file) {
  con <- file(gmt_file, open = "wt")
  on.exit(close(con), add = TRUE)

  for (i in seq_along(gene_set)) {
    cat(names(gene_set)[i], "\tNA\t", paste(gene_set[[i]], collapse = "\t"), "\n",
        file = con, sep = "")
  }
}

write_gmt(io_gene_list, io_gmt_file)

# =========================================================
# 3. Prepare TCGA LUAD TPM matrix
# =========================================================
expr_data <- fread(tpm_file, header = TRUE, sep = "\t")
expr_data <- as.data.frame(expr_data)

gene.code.v22 <- read.csv(gene_anno_file, header = TRUE, sep = ",")

expr_data[, 1] <- gene.code.v22[match(expr_data[, 1], gene.code.v22$X), ][, 2]
expr_data <- expr_data[!duplicated(expr_data[, 1]), ]
colnames(expr_data) <- substring(colnames(expr_data), 1, 12)
rownames(expr_data) <- expr_data[, 1]
expr_data <- expr_data[, -1, drop = FALSE]
colnames(expr_data) <- gsub(".", "-", colnames(expr_data), fixed = TRUE)

protein.coding.genes.id <- gene.code.v22[gene.code.v22$gene_type == "protein_coding", ]$gene_name
expr_data <- expr_data[rownames(expr_data) %in% protein.coding.genes.id, , drop = FALSE]

fwrite(
  expr_data,
  file = file.path(raw_dir, "tcga_luad_mRNA_TPM.csv"),
  row.names = TRUE,
  quote = TRUE,
  sep = ","
)

# =========================================================
# 4. Read gene sets and run ssGSEA
# =========================================================
io.can.cell.gene_set <- read.gmt(io_gmt_file)
io.can.cell.gene_set$term <- paste0("IO_", io.can.cell.gene_set$term)

hallmark.gene_set <- read.gmt(hallmark_gmt)
kegg.gene_set     <- read.gmt(kegg_gmt)
reactome.gene_set <- read.gmt(reactome_gmt)
biocarta.gene_set <- read.gmt(biocarta_gmt)

comb.gene.set <- rbind.data.frame(
  io.can.cell.gene_set,
  hallmark.gene_set,
  kegg.gene_set,
  reactome.gene_set,
  biocarta.gene_set
)

gene.list <- split(as.matrix(comb.gene.set)[, 2], comb.gene.set[, 1])

tcga.luad.tpm.rawdata <- as.matrix(expr_data)

GSVA.results <- gsva(
  tcga.luad.tpm.rawdata,
  gene.list,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = TRUE
)

write.csv(
  GSVA.results,
  file = file.path(result_dir, "tcga.luad.gsva.results.csv")
)

# =========================================================
# 5. Load cluster labels and prepare ordered annotation
# =========================================================
GSVA.results <- read.csv(
  file = file.path(result_dir, "tcga.luad.gsva.results.csv"),
  header = TRUE,
  row.names = 1,
  as.is = TRUE,
  check.names = FALSE
)

tcga.luad_sam.inf_Cluster_res <- read.csv(
  file = cluster_file,
  header = TRUE,
  as.is = TRUE
)

inter.sam.id <- intersect(tcga.luad_sam.inf_Cluster_res$X, colnames(GSVA.results))
tcga.luad_sam.inf_Cluster_res <- tcga.luad_sam.inf_Cluster_res[
  match(inter.sam.id, tcga.luad_sam.inf_Cluster_res$X),
]

tcga.luad_sam.inf_Cluster_res$consensusClass <- factor(
  tcga.luad_sam.inf_Cluster_res$consensusClass,
  levels = c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4")
)

tcga.luad_sam.inf_Cluster_res_sorted <- tcga.luad_sam.inf_Cluster_res[
  order(tcga.luad_sam.inf_Cluster_res$consensusClass),
]

annotation_c2 <- as.data.frame(
  tcga.luad_sam.inf_Cluster_res_sorted[, c("consensusClass", "cancer.group")]
)
rownames(annotation_c2) <- tcga.luad_sam.inf_Cluster_res_sorted$X

ann_colors2 <- list(
  consensusClass = c(
    "Cluster_1" = "#FFF200",
    "Cluster_2" = "#C0504D",
    "Cluster_3" = "#4F81BD",
    "Cluster_4" = "#9BBB59"
  ),
  cancer.group = c("LUAD" = "red", "LUSC" = "blue"),
  sig.type = c("HALLMARK" = "#1F78B4", "IO" = "#B2DF8A", "KEGG" = "#A6CEE3")
)

split_names <- data.frame(
  Prefix = sub("_.*", "", rownames(GSVA.results)),
  Suffix = sub("^[^_]*_", "", rownames(GSVA.results)),
  stringsAsFactors = FALSE
)

annotation_row <- cbind.data.frame(sig.type = split_names$Prefix)
rownames(annotation_row) <- rownames(GSVA.results)

# =========================================================
# 6. Plot GSVA heatmap
# =========================================================
p.data <- as.matrix(GSVA.results[, tcga.luad_sam.inf_Cluster_res_sorted$X])

heatmap_obj <- pheatmap(
  p.data,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_method = "mcquitty",
  color = c(
    colorRampPalette(colors = c("#3F48CC", "#809EF7"))(40),
    colorRampPalette(colors = c("#809EF7", "white"))(25),
    colorRampPalette(colors = c("white", "#D55043"))(20),
    colorRampPalette(colors = c("#D55043", "#560b11"))(40)
  ),
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize = 1,
  display_numbers = FALSE,
  annotation_col = annotation_c2,
  annotation_colors = ann_colors2,
  annotation_row = annotation_row
)

pdf(file.path(plot_dir, "tcga.luad.4.cluster.heatmap.pdf"), width = 12, height = 5)
print(heatmap_obj)
dev.off()

# =========================================================
# 7. Test GSVA pathway differences across four clusters
# =========================================================
merged_data <- cbind.data.frame(
  tcga.luad_sam.inf_Cluster_res[match(inter.sam.id, tcga.luad_sam.inf_Cluster_res$X), ][, 1:2],
  t(GSVA.results[, inter.sam.id])
)
rownames(merged_data) <- merged_data$X
colnames(merged_data)[1] <- "ID"

merged_data$consensusClass <- factor(
  merged_data$consensusClass,
  levels = c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4")
)

results <- data.frame()
pathways <- colnames(merged_data)[-(1:2)]

for (pathway in pathways) {
  for (cluster in levels(merged_data$consensusClass)) {
    data_cluster <- merged_data[merged_data$consensusClass == cluster, pathway]
    data_others  <- merged_data[merged_data$consensusClass != cluster, pathway]

    mean_cluster <- mean(data_cluster, na.rm = TRUE)
    mean_others  <- mean(data_others, na.rm = TRUE)
    mean_fold_change <- mean_cluster / mean_others

    test_result <- t.test(data_cluster, data_others, var.equal = FALSE)

    results <- rbind(
      results,
      data.frame(
        Pathway = pathway,
        Cluster = cluster,
        MeanFoldChange = mean_fold_change,
        Tvalue = unname(test_result$statistic),
        Pvalue = test_result$p.value
      )
    )
  }
}

results$fdr <- p.adjust(results$Pvalue, method = "fdr")

diff.res.comb <- cbind.data.frame(
  data.frame(
    Set = sub("_.*", "", results$Pathway),
    Term = sub("^[^_]*_", "", results$Pathway),
    stringsAsFactors = FALSE
  ),
  results
)

write.csv(
  diff.res.comb,
  file = file.path(result_dir, "tcga.luad.gsva.4group.t.test.diff.results.csv"),
  row.names = FALSE
)

network.res <- diff.res.comb[diff.res.comb$fdr < 0.05, ]
write.csv(
  network.res,
  file = file.path(cyto_dir, "network.res.csv"),
  row.names = FALSE
)

results_count <- diff.res.comb %>%
  filter(fdr < 0.05) %>%
  group_by(Term, Set, Pathway) %>%
  summarise(Count = n(), .groups = "drop")

cluster_num <- cbind.data.frame(
  Term = c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4"),
  Set = c("Cluster", "Cluster", "Cluster", "Cluster"),
  Pathway = c(10, 10, 10, 10),
  Count = c(10, 10, 10, 10)
)

node.res <- rbind.data.frame(results_count, cluster_num)
write.csv(
  node.res,
  file = file.path(cyto_dir, "node.res.csv"),
  row.names = FALSE
)

# =========================================================
# 8. Compare CN-derived clusters with published IO subtype
# =========================================================
IO_subtype.inf <- read.csv(io_subtype_inf_file, header = TRUE, as.is = TRUE)
LUAD.IO_subtype.inf <- IO_subtype.inf[IO_subtype.inf$TCGA_project == "LUAD", ]

tcga.sam.inf.comb <- merge(
  tcga.luad_sam.inf_Cluster_res,
  LUAD.IO_subtype.inf,
  by.x = "X",
  by.y = "X",
  all.x = TRUE
)

table_data <- table(tcga.sam.inf.comb$consensusClass, tcga.sam.inf.comb$MFP)
chi_sq_test <- chisq.test(table_data)
print(chi_sq_test)

tcga.sam.inf.comb.new <- tcga.sam.inf.comb[!is.na(tcga.sam.inf.comb$MFP), ]
print(adjustedRandIndex(tcga.sam.inf.comb.new$consensusClass, tcga.sam.inf.comb.new$MFP))
print(CramerV(table_data))

allu.sam.inf <- tcga.sam.inf.comb[, c("X", "consensusClass", "MFP")]
colnames(allu.sam.inf)[1] <- "sample.id"

LUAD.Data <- group_by(allu.sam.inf, consensusClass, MFP) %>%
  summarise(count = n(), .groups = "drop")
LUAD.Data <- cbind.data.frame(cluster.group.1 = LUAD.Data$consensusClass, LUAD.Data)

LUAD_long <- to_lodes_form(
  data.frame(LUAD.Data),
  key = "Demographic",
  axes = 2:3
)
LUAD_long$Demographic <- factor(LUAD_long$Demographic, levels = c("consensusClass", "MFP"))

sankey.plot <- ggplot(
  data = LUAD_long,
  aes(
    x = Demographic,
    stratum = stratum,
    alluvium = alluvium,
    y = count,
    label = stratum
  )
) +
  geom_alluvium(aes(fill = cluster.group.1), alpha = 0.5) +
  scale_fill_manual(values = c(
    "Cluster_1" = "#FFF200",
    "Cluster_2" = "#C0504D",
    "Cluster_3" = "#4F81BC",
    "Cluster_4" = "#9BBB59"
  )) +
  geom_stratum() +
  geom_text(stat = "stratum") +
  theme_classic()

pdf(file.path(sankey_dir, "tcga.luad_subtypes.comp_sankey.plot.pdf"), width = 3.97, height = 10.41)
print(sankey.plot)
dev.off()

# =========================================================
# 9. Cluster-specific GSEA using limma ranking
# =========================================================
tcga_luad_mRNA_TPM <- read.csv(
  file = file.path(raw_dir, "tcga_luad_mRNA_TPM.csv"),
  header = TRUE,
  row.names = 1,
  as.is = TRUE,
  check.names = FALSE
)

inter.sam.id <- intersect(colnames(tcga_luad_mRNA_TPM), tcga.luad_sam.inf_Cluster_res$X)
tcga_luad_mRNA_TPM.new <- tcga_luad_mRNA_TPM[, inter.sam.id, drop = FALSE]
tcga.luad_sam.inf_Cluster_res.new <- tcga.luad_sam.inf_Cluster_res[
  match(inter.sam.id, tcga.luad_sam.inf_Cluster_res$X),
]

comb.gene.site <- rbind(
  io.can.cell.gene_set,
  hallmark.gene_set,
  kegg.gene_set,
  reactome.gene_set,
  biocarta.gene_set
)

for (i in 1:4) {
  cluster_ids <- tcga.luad_sam.inf_Cluster_res.new[
    tcga.luad_sam.inf_Cluster_res.new$consensusClass == paste0("Cluster_", i),
  ]$X

  other_ids <- setdiff(tcga.luad_sam.inf_Cluster_res.new$X, cluster_ids)
  all_ids <- c(other_ids, cluster_ids)

  group_list <- factor(c(rep("nt", length(other_ids)), rep("t", length(cluster_ids))))
  design <- model.matrix(~ 1 + group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- all_ids

  expr_sub <- log2(tcga_luad_mRNA_TPM.new + 1)[, all_ids, drop = FALSE]

  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit, trend = TRUE)
  output <- topTable(fit, coef = 2, n = Inf)

  pre.rank.res <- output
  pre.rank.res$fcSign <- sign(pre.rank.res$logFC)
  pre.rank.res$logP <- -log10(pmax(pre.rank.res$P.Value, .Machine$double.xmin))
  pre.rank.res$metric <- pre.rank.res$logP * pre.rank.res$fcSign

  y <- cbind.data.frame(genename = rownames(pre.rank.res), metric = pre.rank.res$metric)
  z <- y[complete.cases(y), ]
  geneList <- setNames(sort(z$metric, decreasing = TRUE), z$genename)

  gsea.res <- GSEA(
    geneList,
    TERM2GENE = comb.gene.site,
    pvalueCutoff = 1,
    nPerm = 1000,
    minGSSize = 1,
    pAdjustMethod = "BH",
    seed = 1
  )

  gsea.res.1 <- gsea.res@result
  write.csv(gsea.res.1, file = file.path(gsea_dir, paste0("cluster", i, ".gsea.res.csv")))
}

# =========================================================
# 10. Combine significant positive GSEA results
# =========================================================
cluster1.gsea.res <- cbind.data.frame(
  Cluster = "cluster1",
  read.csv(file.path(gsea_dir, "cluster1.gsea.res.csv"), header = TRUE, row.names = 1, as.is = TRUE)
)
cluster2.gsea.res <- cbind.data.frame(
  Cluster = "cluster2",
  read.csv(file.path(gsea_dir, "cluster2.gsea.res.csv"), header = TRUE, row.names = 1, as.is = TRUE)
)
cluster3.gsea.res <- cbind.data.frame(
  Cluster = "cluster3",
  read.csv(file.path(gsea_dir, "cluster3.gsea.res.csv"), header = TRUE, row.names = 1, as.is = TRUE)
)
cluster4.gsea.res <- cbind.data.frame(
  Cluster = "cluster4",
  read.csv(file.path(gsea_dir, "cluster4.gsea.res.csv"), header = TRUE, row.names = 1, as.is = TRUE)
)

cluster1.gsea.res.sig <- cluster1.gsea.res[cluster1.gsea.res$pvalue < 0.1 & cluster1.gsea.res$NES > 0, ][1:300, ]
cluster2.gsea.res.sig <- cluster2.gsea.res[cluster2.gsea.res$pvalue < 0.1 & cluster2.gsea.res$NES > 0, ][1:300, ]
cluster3.gsea.res.sig <- cluster3.gsea.res[cluster3.gsea.res$pvalue < 0.1 & cluster3.gsea.res$NES > 0, ][1:300, ]
cluster4.gsea.res.sig <- cluster4.gsea.res[cluster4.gsea.res$pvalue < 0.1 & cluster4.gsea.res$NES > 0, ][1:300, ]

cluster.gsea.sig.comb.res <- rbind.data.frame(
  cluster1.gsea.res.sig,
  cluster2.gsea.res.sig,
  cluster3.gsea.res.sig,
  cluster4.gsea.res.sig
)

cluster.gsea.sig.comb.res <- cbind.data.frame(
  data.frame(
    Set = sub("_.*", "", cluster.gsea.sig.comb.res$ID),
    Term = sub("^[^_]*_", "", cluster.gsea.sig.comb.res$ID),
    stringsAsFactors = FALSE
  ),
  cluster.gsea.sig.comb.res
)

write.csv(
  cluster.gsea.sig.comb.res,
  file = file.path(gsea_dir, "cluster.gsea.sig.comb.res.csv"),
  row.names = FALSE
)

# =========================================================
# 11. Cluster-specific highlighted bubble plots
# =========================================================
make_cluster_bubble <- function(df, cluster_id, labels_keep, fill_color, reverse_x = FALSE, reverse_y = FALSE) {
  each.cluster.gsea.sig.res <- df[df$Cluster == cluster_id, ]
  each.cluster.gsea.sig.res <- each.cluster.gsea.sig.res %>%
    arrange(pvalue, desc(NES))

  other_ids <- setdiff(each.cluster.gsea.sig.res$ID, labels_keep)
  other_ids <- head(other_ids, 50 - length(labels_keep))

  gsea.data <- rbind.data.frame(
    each.cluster.gsea.sig.res[, c("Set", "ID", "NES", "setSize", "pvalue")][match(other_ids, each.cluster.gsea.sig.res$ID), ],
    each.cluster.gsea.sig.res[, c("Set", "ID", "NES", "setSize", "pvalue")][match(labels_keep, each.cluster.gsea.sig.res$ID), ]
  )

  gsea.data$setSize_1 <- gsea.data$setSize / 10
  gsea.data <- gsea.data[order(gsea.data$NES, decreasing = TRUE), ]
  gsea.data$ID <- factor(gsea.data$ID, levels = gsea.data$ID)
  gsea.data$xlab <- seq_len(nrow(gsea.data))

  data_label <- gsea.data[gsea.data$ID %in% labels_keep, ]

  p <- ggplot(data = gsea.data, aes(x = xlab, y = NES)) +
    geom_point(
      aes(size = setSize_1, alpha = -log10(pmax(pvalue, .Machine$double.xmin))),
      shape = 21,
      stroke = 0.7,
      fill = fill_color,
      colour = "black"
    ) +
    scale_size_continuous(range = c(1, 6)) +
    xlab("Gene sets") +
    ylab("Normalized enrichment score (NES)") +
    theme_classic(base_size = 15) +
    guides(
      size = guide_legend(title = "Detection\n(proportion)"),
      alpha = guide_legend(title = "Significance\n(-log10 p-val.)")
    ) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.text = element_text(face = "bold"),
      axis.title = element_text(size = 5)
    )

  if (reverse_x) p <- p + scale_x_reverse()
  if (reverse_y) p <- p + scale_y_reverse()

  p +
    geom_text_repel(
      data = data_label,
      aes(label = ID),
      color = "black",
      size = 1,
      force = 10,
      point.padding = 0.5,
      min.segment.length = 0,
      segment.color = "grey20",
      segment.size = 0.3,
      segment.alpha = 0.8,
      nudge_y = -0.1
    )
}

cluster1.anno <- c(
  "REACTOME_TICAM1_TRAF6_DEPENDENT_INDUCTION_OF_TAK1_COMPLEX",
  "BIOCARTA_TOLL_PATHWAY",
  "REACTOME_NUCLEAR_ENVELOPE_NE_REASSEMBLY",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
  "REACTOME_RHO_GTPASES_ACTIVATE_FORMINS",
  "BIOCARTA_NFKB_PATHWAY",
  "IO_Tumor proliferation rate",
  "REACTOME_CELL_CYCLE_CHECKPOINTS"
)

cluster2.anno <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "REACTOME_MAPK6_MAPK4_SIGNALING",
  "REACTOME_REGULATION_OF_RUNX2_EXPRESSION_AND_ACTIVITY",
  "REACTOME_NEGATIVE_REGULATION_OF_NOTCH4_SIGNALING"
)

cluster3.anno <- c(
  "REACTOME_CELL_CYCLE_CHECKPOINTS",
  "HALLMARK_G2M_CHECKPOINT",
  "REACTOME_MITOTIC_SPINDLE_CHECKPOINT",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

cluster4.anno <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"
)

cluster1.bubble.plot <- make_cluster_bubble(cluster.gsea.sig.comb.res, "cluster1", cluster1.anno, "#FFF200")
cluster2.bubble.plot <- make_cluster_bubble(cluster.gsea.sig.comb.res, "cluster2", cluster2.anno, "#C0504D", reverse_x = TRUE)
cluster3.bubble.plot <- make_cluster_bubble(cluster.gsea.sig.comb.res, "cluster3", cluster3.anno, "#4F81BD", reverse_y = TRUE)
cluster4.bubble.plot <- make_cluster_bubble(cluster.gsea.sig.comb.res, "cluster4", cluster4.anno, "#9BBB59")

pdf(file.path(gsea_plot_dir, "cluster1.bubble.plot.pdf"), width = 6.81, height = 4.39)
print(cluster1.bubble.plot)
dev.off()

pdf(file.path(gsea_plot_dir, "cluster2.bubble.plot.pdf"), width = 6.81, height = 4.39)
print(cluster2.bubble.plot)
dev.off()

pdf(file.path(gsea_plot_dir, "cluster3.bubble.plot.pdf"), width = 6.81, height = 4.39)
print(cluster3.bubble.plot)
dev.off()

pdf(file.path(gsea_plot_dir, "cluster4.bubble.plot.pdf"), width = 6.81, height = 4.39)
print(cluster4.bubble.plot)
dev.off()

# =========================================================
# 12. Combined cluster GSEA bubble plot
# =========================================================
cluster1.gsea.data.1 <- cluster.gsea.sig.comb.res[cluster.gsea.sig.comb.res$Cluster == "cluster1", ]
cluster2.gsea.data.1 <- cluster.gsea.sig.comb.res[cluster.gsea.sig.comb.res$Cluster == "cluster2", ]
cluster3.gsea.data.1 <- cluster.gsea.sig.comb.res[cluster.gsea.sig.comb.res$Cluster == "cluster3", ]
cluster4.gsea.data.1 <- cluster.gsea.sig.comb.res[cluster.gsea.sig.comb.res$Cluster == "cluster4", ]

cluster1.data_label <- cluster1.gsea.data.1[cluster1.gsea.data.1$ID %in% cluster1.anno, ]
cluster2.data_label <- cluster2.gsea.data.1[cluster2.gsea.data.1$ID %in% cluster2.anno, ]
cluster3.data_label <- cluster3.gsea.data.1[cluster3.gsea.data.1$ID %in% cluster3.anno, ]
cluster4.data_label <- cluster4.gsea.data.1[cluster4.gsea.data.1$ID %in% cluster4.anno, ]

gsea.data.comb <- rbind.data.frame(
  cluster1.gsea.data.1,
  cluster2.gsea.data.1,
  cluster3.gsea.data.1,
  cluster4.gsea.data.1
)

gsea.data.comb$xlab <- ave(gsea.data.comb$NES, gsea.data.comb$Cluster, FUN = seq_along)
gsea.data.comb$clust.and.ID <- paste0(gsea.data.comb$ID, "_", gsea.data.comb$Cluster)

label.comb <- rbind.data.frame(
  cluster1.data_label,
  cluster2.data_label,
  cluster3.data_label,
  cluster4.data_label
)
label.comb$clust.and.ID <- paste0(label.comb$ID, "_", label.comb$Cluster)

gsea.data.comb$alpha_group <- ifelse(
  gsea.data.comb$clust.and.ID %in% label.comb$clust.and.ID,
  "p", "np"
)
gsea.data.comb$log.pvalue <- -log10(pmax(gsea.data.comb$pvalue, .Machine$double.xmin))

each.cluster_gsea_bubble.plot <- ggplot(data = gsea.data.comb, aes(x = xlab, y = NES)) +
  geom_point(
    aes(size = log.pvalue, alpha = alpha_group, fill = Cluster, color = alpha_group),
    shape = 21
  ) +
  scale_size_continuous(range = c(0.2, 10)) +
  scale_alpha_manual(values = c("np" = 0.1, "p" = 0.9)) +
  scale_color_manual(values = c("np" = NA, "p" = "black")) +
  scale_fill_manual(values = c(
    "cluster1" = "#EFBC20",
    "cluster2" = "#C0504D",
    "cluster3" = "#4F81BC",
    "cluster4" = "#9BBB59"
  )) +
  theme_classic() +
  guides(
    size = guide_legend(title = "Significance\n(-log10 p-val.)"),
    alpha = guide_legend(title = "Detection"),
    fill = guide_legend(title = "Cluster")
  ) +
  facet_wrap(~ Cluster, scales = "free_x") +
  geom_text_repel(
    data = label.comb,
    aes(label = ID, color = Cluster),
    size = 3,
    force = 10,
    point.padding = 0.5,
    segment.color = "grey20",
    segment.size = 0.3,
    segment.alpha = 0.8,
    nudge_y = -0.1
  )

pdf(file.path(gsea_plot_dir, "each.cluster_gsea_bubble.plot.pdf"), width = 8.78, height = 10.95)
print(each.cluster_gsea_bubble.plot)
dev.off()

write.csv(gsea.data.comb, file = file.path(gsea_dir, "cluster_4group_gsea.data.comb.csv"), row.names = FALSE)

# =========================================================
# 13. External IO-related quantitative signatures
# =========================================================
TCGA_pan.cancer_IO_sig.score <- read.csv(
  file = io_quant_score_file,
  header = TRUE,
  row.names = 1,
  as.is = TRUE
)

LUAD.io.sig.score <- TCGA_pan.cancer_IO_sig.score[
  TCGA_pan.cancer_IO_sig.score$TCGA.Study == "LUAD",
]

inter.sam.id <- intersect(tcga.luad_sam.inf_Cluster_res$X, TCGA_pan.cancer_IO_sig.score$Sample.ID)

tcga.luad_cluster_io.sig <- cbind.data.frame(
  tcga.luad_sam.inf_Cluster_res[match(inter.sam.id, tcga.luad_sam.inf_Cluster_res$X), ],
  TCGA_pan.cancer_IO_sig.score[match(inter.sam.id, TCGA_pan.cancer_IO_sig.score$Sample.ID), ]
)

tcga.luad_cluster_io.sig_selection <- tcga.luad_cluster_io.sig[, c(
  "consensusClass", "age_at_initial_pathologic_diagnosis", "gender",
  "Intratumor.Heterogeneity", "SNV.Neoantigens", "Aneuploidy.Score",
  "TMB", "MSI", "INFR.score", "EIGS.score", "TLS.score",
  "Stroma_EMT_TGFb", "Angiogenesis"
)]
rownames(tcga.luad_cluster_io.sig_selection) <- tcga.luad_cluster_io.sig$X
colnames(tcga.luad_cluster_io.sig_selection)[2] <- "age"

s.io.id <- c(
  "consensusClass", "Intratumor.Heterogeneity", "SNV.Neoantigens",
  "Aneuploidy.Score", "TMB", "MSI", "INFR.score", "EIGS.score",
  "TLS.score", "Stroma_EMT_TGFb", "Angiogenesis"
)

g.data.raw <- tcga.luad_cluster_io.sig_selection[, s.io.id]
g.data.raw.1 <- reshape2::melt(g.data.raw, value.name = "consensusClass")
colnames(g.data.raw.1)[3] <- "value"

io_violin_plot <- ggplot(g.data.raw.1, aes(x = consensusClass, y = value)) +
  geom_violin(
    aes(fill = consensusClass),
    scale = "area",
    alpha = 0.8,
    trim = FALSE,
    adjust = 1,
    color = "black"
  ) +
  facet_wrap(~ variable, scales = "free_y", nrow = 5, ncol = 2, strip.position = "top") +
  labs(x = "Consensus Class", y = "Signature value") +
  theme_classic() +
  scale_fill_manual(values = c(
    "Cluster_1" = "#EFBC20",
    "Cluster_2" = "#C0504D",
    "Cluster_3" = "#4F81BC",
    "Cluster_4" = "#9BBB59"
  ))

data_summary <- function(x) {
  m <- mean(x, na.rm = TRUE)
  ymin <- m - sd(x, na.rm = TRUE)
  ymax <- m + sd(x, na.rm = TRUE)
  c(y = m, ymin = ymin, ymax = ymax)
}

io_violin_plot2 <- io_violin_plot +
  stat_summary(
    fun.data = data_summary,
    geom = "pointrange",
    color = "black",
    linewidth = 0.5
  )

pdf(file.path(gsea_plot_dir, "each.cluster_IO.sig_plot.pdf"), width = 6.34, height = 8.05)
print(io_violin_plot2)
dev.off()

# =========================================================
# 14. ANOVA for external IO signatures across clusters
# =========================================================
tcga.luad_cluster_io.sig_selection$consensusClass <- factor(tcga.luad_cluster_io.sig_selection$consensusClass)
tcga.luad_cluster_io.sig_selection$gender <- factor(tcga.luad_cluster_io.sig_selection$gender)

s.io.id <- c(
  "Intratumor.Heterogeneity", "SNV.Neoantigens", "Aneuploidy.Score",
  "TMB", "MSI", "INFR.score", "EIGS.score", "TLS.score",
  "Stroma_EMT_TGFb", "Angiogenesis"
)

io_sig_results <- data.frame(Pathway = character(), Pvalue = numeric(), Fvalue = numeric(), stringsAsFactors = FALSE)

for (pathway in s.io.id) {
  formula <- as.formula(paste(pathway, "~ consensusClass + age + gender"))
  model <- lm(formula, data = tcga.luad_cluster_io.sig_selection)

  anova_results <- Anova(model, type = "II")
  anova_consensusClass <- anova_results["consensusClass", ]

  io_sig_results <- rbind(
    io_sig_results,
    data.frame(
      Pathway = pathway,
      Pvalue = anova_consensusClass$`Pr(>F)`,
      Fvalue = anova_consensusClass$`F value`
    )
  )
}

io_sig_results$Significance <- ifelse(
  io_sig_results$Pvalue < 0.001, "***",
  ifelse(io_sig_results$Pvalue < 0.01, "**",
         ifelse(io_sig_results$Pvalue < 0.05, "*", "n.s."))
)

write.csv(
  io_sig_results,
  file = file.path(gsea_dir, "IO_sig_4group_diff.res.csv"),
  row.names = FALSE
)