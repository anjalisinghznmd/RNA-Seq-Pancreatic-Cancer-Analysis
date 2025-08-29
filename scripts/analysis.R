# -------------------------------------------------------------------
# RNA-seq Differential Gene Expression Analysis
# Project: Pancreatic Cancer Metastasis (GSE158527)
# FINAL CORRECTED VERSION (with built-in metadata)
# -------------------------------------------------------------------

# --- 1. Load Required Libraries ---
library(DESeq2)
library(ggplot2)
library(pheatmap)

# --- 2. Load the Counts Data ---
counts_data <- read.csv("data/GSE158527_raw_counts.csv", row.names = 1)

# --- 3. Create Metadata Directly in R ---
# This step bypasses any issues with reading the metadata.csv file.
# We create the table manually here. Note the use of dots (.) instead of hyphens (-).
sample_ids <- c("S2.013.1", "S2.013.2", "S2.013.3",
                "S2.LM7AA.1", "S2.LM7AA.2", "S2.LM7AA.3",
                "Suit2.1", "Suit2.2", "Suit2.3",
                "Suit2.OE.1", "Suit2.OE.2", "Suit2.OE.3")

conditions <- c("S2-013", "S2-013", "S2-013",
                "Metastatic", "Metastatic", "Metastatic",
                "Parental", "Parental", "Parental",
                "Suit2-OE", "Suit2-OE", "Suit2-OE")

# Create the metadata data frame.
metadata <- data.frame(Condition = conditions, row.names = sample_ids)


# --- 4. Data Preparation and Sanity Checks ---
# Let's ensure the metadata is in the same order as the counts data.
metadata <- metadata[colnames(counts_data), , drop = FALSE]

# This check MUST return TRUE now.
cat("Checking if column names and row names match...\n")
print(all(colnames(counts_data) == rownames(metadata)))


# --- 5. Create the DESeq2 Dataset Object ---
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = metadata,
                              design = ~ Condition)

# --- 6. Pre-filtering ---
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# --- 7. Set the Reference Level ---
dds$Condition <- relevel(dds$Condition, ref = "Parental")

# --- 8. Run the DESeq2 Analysis ---
cat("Running DESeq2 analysis... This may take a minute.\n")
dds <- DESeq(dds)

# --- 9. Extract the Results ---
cat("Extracting results...\n")
res <- results(dds, contrast = c("Condition", "Metastatic", "Parental"))
res_df <- as.data.frame(res)
res_df <- na.omit(res_df)
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, "results/pancreatic_cancer_DE_results.csv")
cat("Top 10 Differentially Expressed Genes:\n")
print(head(res_df, 10))

# --- 10. Generate a Volcano Plot ---
cat("Generating volcano plot...\n")
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot: Metastatic vs. Parental Cells",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-value)") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40")

ggsave("results/volcano_plot.png", plot = volcano_plot, width = 10, height = 8)
cat("Volcano plot saved to results/volcano_plot.png\n")

# --- 11. Generate a Heatmap (FINAL Corrected Version) ---
cat("Generating heatmap...\n")

# Get the names of the top 50 most significant genes from our results.
top_50_genes <- rownames(head(res_df, 50))

# Normalize the counts data for visualization purposes (VST).
vst_counts <- vst(dds, blind = FALSE)

# Robustly select the names of the samples we are interested in.
sample_names <- rownames(colData(dds)[colData(dds)$Condition %in% c("Metastatic", "Parental"), ])

# Subset the transformed counts to only our top 50 genes and selected samples.
heatmap_data <- assay(vst_counts)[top_50_genes, sample_names]

# Create the annotation data frame for pheatmap. This is the key fix.
# It needs to be a data frame with row names that match the columns of heatmap_data.
annotation_df <- as.data.frame(colData(dds)[sample_names, "Condition", drop = FALSE])
colnames(annotation_df) <- "Condition" # Ensure the column is named correctly

# Create the heatmap and save it.
png("results/heatmap.png", width = 800, height = 1000)
pheatmap(heatmap_data,
         scale = "row",
         annotation_col = annotation_df,
         main = "Heatmap of Top 50 DE Genes",
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = TRUE)
dev.off()

cat("Heatmap saved to results/heatmap.png\n")
cat("Analysis complete!\n")
