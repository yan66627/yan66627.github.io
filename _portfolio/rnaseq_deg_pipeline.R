# =============================================================================
# RNA-seq Differential Expression Analysis Pipeline
# Author: Yan Jin
#
# Merges two GEO count matrices, applies CPM-based filtering, runs DESeq2
# differential expression (DHT vs ctrl), and produces a volcano plot with
# labeled up/down-regulated genes.
#
# Data: GSE220168 + GSE205885 (merged AR target gene RNA-seq)
# Comparison: 16h DHT treatment vs control
# Dependencies: DESeq2, ggplot2, RColorBrewer
# =============================================================================

library(DESeq2)
library(ggplot2)
library(RColorBrewer)

# -----------------------------------------------------------------------------
# 1. Load and merge count matrices from two GEO datasets
# -----------------------------------------------------------------------------
load_counts <- function(path) {
  n_cols <- length(read.table(path, sep = "\t", nrows = 1))
  read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
             colClasses = c("character", rep("numeric", n_cols - 1)))
}

cts1 <- load_counts("GSE220168_merged_gene_name_expression.txt")
cts2 <- load_counts("GSE205885_merged_gene_name_expression.txt")

gene_names <- cts1[, 1]
cts <- cbind(cts1[, -c(1:2)], cts2[, -c(1:2)])
rownames(cts) <- gene_names

# -----------------------------------------------------------------------------
# 2. Load sample metadata and align with count matrix
# -----------------------------------------------------------------------------
coldata <- read.csv("coldata_rna_seq_filtered.csv", row.names = "X")
rownames(coldata) <- coldata$SRR
coldata$SRR <- NULL

# Reorder columns to match metadata row order
cts <- cts[, rownames(coldata), drop = FALSE]
stopifnot(all(rownames(coldata) == colnames(cts)))

# -----------------------------------------------------------------------------
# 3. Build DESeq2 object with CPM-based pre-filtering
#    Keep genes with CPM > 1 in at least 3 samples
# -----------------------------------------------------------------------------
design_mat <- model.matrix(~ 0 + factor(coldata$Time_treatment))
colnames(design_mat) <- levels(factor(coldata$Time_treatment))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData   = coldata,
                              design    = design_mat)

cpm    <- apply(cts, 2, function(x) x / sum(x) * 1e6)
keep   <- rowSums(cpm > 1) >= 3
dds    <- dds[keep, ]
cat("Genes retained after CPM filtering:", sum(keep), "\n")

# -----------------------------------------------------------------------------
# 4. Run DESeq2 and extract results (16h DHT vs ctrl)
# -----------------------------------------------------------------------------
dds <- DESeq(dds)
res <- results(dds, contrast = list("X16h_DHT", "ctrl"))
res$padj[is.na(res$padj)] <- 1

res_df <- as.data.frame(res)

# Filter: FDR < 0.05 and |log2FC| > 1
res_filtered <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
res_up        <- res_filtered[res_filtered$log2FoldChange >  1, ]
res_down      <- res_filtered[res_filtered$log2FoldChange < -1, ]

cat("Total DEGs (padj < 0.05 & |log2FC| > 1):", nrow(res_filtered), "\n")
cat("Upregulated:", nrow(res_up), "| Downregulated:", nrow(res_down), "\n")

# -----------------------------------------------------------------------------
# 5. Save results
# -----------------------------------------------------------------------------
write.csv(res_df,       "DHT_vs_ctrl_all_results.csv")
write.csv(res_filtered, "DHT_vs_ctrl_filtered_DEG.csv")
write.csv(res_up,       "DHT_vs_ctrl_UP.csv")
write.csv(res_down,     "DHT_vs_ctrl_DOWN.csv")

# -----------------------------------------------------------------------------
# 6. Volcano plot
# -----------------------------------------------------------------------------
res_df$significance <- "Not Sig"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "Up"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significance), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title  = "Volcano Plot: 16h DHT vs ctrl",
       x      = "log2 Fold Change",
       y      = "-log10(p-value)",
       color  = "Significance") +
  theme_minimal(base_size = 14)
