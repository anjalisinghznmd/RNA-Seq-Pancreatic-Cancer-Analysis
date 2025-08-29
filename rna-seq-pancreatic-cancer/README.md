RNA-seq Analysis of Pancreatic Cancer Metastasis
Project Overview
This project performs a comprehensive differential gene expression analysis on a public RNA-sequencing (RNA-seq) dataset to identify genes associated with metastasis in pancreatic cancer. The primary goal is to compare the gene expression profiles of metastatic cancer cells against their less aggressive, parental counterparts to uncover the molecular drivers of cancer progression.

This is a foundational bioinformatics project that demonstrates core skills in data curation, statistical analysis, and data visualization using industry-standard tools like R and Bioconductor.

Dataset
The dataset used in this analysis is publicly available from the NCBI Gene Expression Omnibus (GEO) database.

GEO Accession: GSE158527

Title: Glycosyltransferase ST6Gal-I promotes the epithelial to mesenchymal transition in pancreatic cancer cells

Organism: Homo sapiens (Human)

Experimental Groups: The analysis focuses on comparing the "Suit2 parental pancreatic cancer cells" (Parental) against the "Suit2-LM7AA metastatic pancreatic cancer cells" (Metastatic).

Analysis Workflow & Tools
The analysis was conducted entirely in R, leveraging powerful packages from the Bioconductor project.

Data Curation & Preparation:

Raw gene counts were downloaded from the GEO supplementary files.

A metadata table was created to map sample IDs to their respective experimental conditions (Parental vs. Metastatic).

Differential Expression Analysis:

The DESeq2 package was used to perform the statistical analysis.

The workflow included normalization of raw counts, estimation of data dispersion, and fitting a negative binomial model to identify genes with statistically significant changes in expression.

Visualization:

ggplot2 was used to create a volcano plot to visualize the relationship between fold change and statistical significance for all genes.

pheatmap was used to generate a heatmap of the top 50 differentially expressed genes, revealing distinct expression patterns between the two experimental groups.

Core Tools:

R: The primary programming language for the analysis.

DESeq2: The statistical engine for determining differential expression.

ggplot2 & pheatmap: For generating publication-quality visualizations.

Key Results
The analysis successfully identified a large number of genes that are differentially expressed between the metastatic and parental cancer cell lines.

1. Volcano Plot
The volcano plot below provides a comprehensive overview of the analysis results. Each dot represents a single gene. The red dots highlight genes that are both statistically significant (Adjusted P-value < 0.05) and have a substantial change in expression (Log2 Fold Change > 1 or < -1).

This plot clearly shows numerous genes that are strongly upregulated (right side) and downregulated (left side) in the metastatic cells, providing a rich set of candidates for further biological investigation.

2. Heatmap of Top 50 Genes
The heatmap visualizes the expression levels of the 50 most statistically significant genes across the samples. The samples are grouped by condition (Metastatic vs. Parental).

This heatmap reveals two distinct gene expression signatures:

A cluster of genes that are highly expressed (red) in the Metastatic samples but have low expression (blue) in the Parental samples (upregulated genes).

A second cluster of genes with the opposite pattern—low expression in Metastatic samples and high expression in Parental samples (downregulated genes).

This clear separation demonstrates that the identified genes create a robust molecular signature that distinguishes the aggressive, metastatic cells from their parental counterparts.

How to Reproduce the Analysis
Ensure you have R and RStudio installed.

Clone this repository to your local machine.

Open the scripts/analysis.R script in RStudio.

Run the script. It will automatically install required packages and perform the full analysis, generating all result files in the results/ directory.