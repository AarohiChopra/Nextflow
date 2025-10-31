# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)
library(argparse)
library(pheatmap)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(stringr)
library(ggplot2)

parser <- ArgumentParser(description = "Run DESeq2 differential expression analysis")

parser$add_argument("--counts", required = TRUE, help = "Path to the counts CSV file")
parser$add_argument("--condition", required = TRUE, help = "Path to the condition TXT file)")
parser$add_argument("--batch", required = TRUE, help = "Path to the batch CSV file")
parser$add_argument("--chemical", required = TRUE, help = "Chemical/treatment name")
parser$add_argument("--conc", required = TRUE, help = "Concentration value")
parser$add_argument("--control", required = TRUE, help = "Control condition label")
parser$add_argument("--outdir", default = ".", help = "Output directory")

args <- parser$parse_args()

safe_label <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
args$conc     <- gsub("\\[|\\]", "", args$conc)  # strip [ ]
args$chemical <- safe_label(args$chemical)
args$conc     <- safe_label(args$conc)
args$control  <- safe_label(args$control)


log_file <- sprintf("DESeq2_log_%s_%s_%s.txt",
                    args$chemical, args$conc, args$control)
sink(file = log_file, split = TRUE)
print("Arguments passed in the run command:")
print(args)

sink(file = log_file, append = TRUE)
counts_data <- read.csv(args$counts)
counts_data <- na.omit(counts_data)
counts_data <- counts_data %>%
  column_to_rownames(var = "Geneid")
if (any(counts_data < 0)) { counts_data[counts_data < 0] <- 0 }
print(tail(counts_data))

condition <- read.table(args$condition, header = TRUE, fill = TRUE, row.names = 1)
batch <- read.csv(args$batch, header = TRUE, fill = TRUE, row.names = 1)

conditions <- as.vector(condition[, 1])
# Check if conditions and counts_data have matching column names
if (ncol(counts_data) != length(conditions)) {
  stop("Error: Number of conditions does not match the number of columns in counts_data")
}
colData <- data.frame(row.names = colnames(counts_data), condition = conditions, batch = batch)  # <-- Added batch to colData
colData$condition <- as.factor(colData$condition)

# Ensure "control" is a valid reference level
#if ("Control" %in% levels(colData$condition)) {
#  colData$condition <- relevel(colData$condition, ref = "Control")
#} else {
#  stop("Error: 'somthing wrong in condition' is not an existing level in the condition column.")
#}

condition_batch <- as.vector(batch[, 1]) # convert to vector
colData$batch <- condition_batch
print("colData: ")
print(colData)
print("Number of rows in the counts data: ")
print(nrow(counts_data))

dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~ condition)
filter2 <- rowSums(counts(dds) >= 10) > (ncol(dds) / 2) # strict filter
dds <- dds[filter2,]
print(dds)

# VST makes variance  approximately constant across the range of mean values.
# Explanation: Reduces the differences in variability (spread) of gene expression levels between genes that are lowly expressed and those that are highly expressed.
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Calculate the Spearman correlation matrix, and keep the samples that meet the cutoff range
cor_matrix <- cor(vsd_mat, method = "spearman")
one_minus_spearman <- 1 - cor_matrix
#print(one_minus_spearman)
cutoff <- 0.1
samples_to_keep <- apply(one_minus_spearman, 2, function(x) any(x < cutoff & x > 0))
samples_names_to_keep <- colnames(one_minus_spearman)[samples_to_keep]
print(samples_names_to_keep)


pdf_file <- paste0("DESeq2_report_", args$chemical, "_", args$conc, "_", args$control, ".pdf")
pdf(pdf_file, onefile = TRUE) 
sample_order <- order(colData$condition)
cor_matrix_copy <- cor_matrix[sample_order, sample_order]
ordered_conditions <- colData$condition[sample_order]

# Replace the row and column names in the correlation matrix with condition names
#rownames(cor_matrix_copy) <- ordered_conditions
#colnames(cor_matrix_copy) <- ordered_conditions
pheatmap(cor_matrix_copy,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         cluster_rows = FALSE,  # Disable clustering
         cluster_cols = FALSE,  # Disable clustering
         color = colorRampPalette(c("white","pink","red","maroon"))(100),
         legend = TRUE,  # Add legend for color
         legend_breaks = seq(-1, 1, 0.2),  # Define breaks for correlation values
         legend_labels = round(seq(-1, 1, 0.2), 2),  # Label the legend
         display_numbers = TRUE,  # Show correlation values on the heatmap
         fontsize_number = 8)  # Adjust font size for the correlation values


sink(file = log_file, append = TRUE)

dds <- DESeq(dds)
res <- results(dds)
summary(res)

significant_genes <- res[which(!is.na(res$padj) & res$padj < 0.05), ]
significant_genes_sorted <- significant_genes[order(significant_genes$padj), ]
print("Significant genes: ")
print(significant_genes_sorted)

# Create a dataframe to create file for gsea
write_res <- results(dds)
df <- data.frame(gene = rownames(write_res), stat = res$log2FoldChange)
df_sorted <- df[order(-df$stat), ]
write.csv(df_sorted, paste0("log2foldsorted_file_", args$chemical, "_", args$conc, "_", args$control, ".csv"), row.names = FALSE)

num_upregulated <- sum(significant_genes$log2FoldChange > 0)
num_downregulated <- sum(significant_genes$log2FoldChange < 0)

cat("Number of upregulated genes:", num_upregulated, "\n")
cat("Number of downregulated genes:", num_downregulated, "\n")

library(ggplot2)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$Treatment <- ifelse(pcaData$condition == "DMSO", "Control", "Treatment")
pcaData$Treatment <- factor(pcaData$Treatment, levels = c("Control", "Treatment"))
pcaData$Concentration <- ifelse(
  pcaData$Treatment == "Control",
  "Control",
  gsub(".*\\[|\\]", "", pcaData$condition)
)
pcaData$Concentration <- as.factor(pcaData$Concentration)
pcaData$SampleID <- colnames(vsd_mat)
pcaData$SampleID <- colnames(vsd)
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Treatment, shape = Concentration, label = SampleID)) +
  geom_point(size = 4) +
  geom_text(vjust = -1.2, size = 3.5, check_overlap = TRUE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(
    values = c("Control" = "#ADD8E6", "Treatment" = "red"),
    name = "Treatment"
  ) +
  scale_shape_manual(
    values = 0:25,
    name = "Concentration"
  ) +
  guides(fill = "none", label = "none", group = "none") +
  theme_minimal(base_size = 16) +
  ggtitle(paste0("PCA: Control vs ", args$chemical,"Treatments"))+
  theme(plot.title = element_text(hjust = 0.5))
options(repr.plot.width = 12, repr.plot.height = 10)
print(p)



# MA plot (log fold change vs mean expression)
plotMA(res, main = "MA Plot", ylim = c(-5, 5))

# PCA plot to visualize sample clustering
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

# Dispersion plot to visualize gene variability
plotDispEsts(dds)

# Cook's distance for outlier detection (for one sample)
plot(seq_along(assays(dds)[["cooks"]][,1]), assays(dds)[["cooks"]][,1], type = "h",
     main = "Cook's Distance for Sample 1", xlab = "Genes", ylab = "Cook's Distance")

# Heatmap of sample distances
sampleDists <- dist(t(vsd_mat))
sampleDistMatrix <- as.matrix(sampleDists)
heatmap(sampleDistMatrix, symm = TRUE, main = "Sample Distance Heatmap")
# the heatmap made above is better so use that

# Close the PDF device
dev.off()



sink(file = log_file, append = TRUE)
normalized_counts <- counts(dds, normalized=TRUE)
valid_rows <- !is.na(rownames(normalized_counts)) & rownames(normalized_counts) != ""
filtered_counts <- normalized_counts[valid_rows, ]
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(filtered_counts),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
filtered_counts <- filtered_counts[!is.na(gene_symbols), ]
rownames(filtered_counts) <- gene_symbols[!is.na(gene_symbols)]
sample_metadata <- as.data.frame(colData(dds))  # ensure plain data.frame
if (!"condition" %in% colnames(sample_metadata)) {
  stop("Metadata does not contain a 'condition' column. Check your sample metadata.")
}

# Keep original labels (e.g., 'DMSO', 'BPA_[50]')
sample_metadata$condition_raw <- as.character(sample_metadata$condition)
sample_metadata$condition_raw <- trimws(sample_metadata$condition_raw)

# Anything exactly equal (case-insensitive) to args$control is control; else treatment
ctrl_label <- tolower(args$control)
is_control <- tolower(sample_metadata$condition_raw) == ctrl_label

sample_metadata$group <- ifelse(is_control, "control", "treatment")

# No NAs expected now
if (any(is.na(sample_metadata$group))) {
  bad <- which(is.na(sample_metadata$group))
  stop(sprintf(
    "Unexpected condition values at samples: %s (condition_raw: %s).",
    paste(rownames(sample_metadata)[bad], collapse = ", "),
    paste(sample_metadata$condition_raw[bad], collapse = ", ")
  ))
}

# Consistent levels (control first)
sample_metadata$group <- factor(sample_metadata$group, levels = c("control","treatment"))

# Order samples control first (tie-break by the original label for readability)
sorted_idx <- order(sample_metadata$group, sample_metadata$condition_raw)
sorted_sample_names <- rownames(sample_metadata)[sorted_idx]

message(
  "Sample order before sorting: ",
  paste(rownames(sample_metadata), as.character(sample_metadata$group), sep=":", collapse=", ")
)
message(
  "Sample order after  sorting: ",
  paste(sorted_sample_names, as.character(sample_metadata$group[sorted_idx]), sep=":", collapse=", ")
)

# Reorder columns in the counts matrix to match
filtered_counts <- filtered_counts[, sorted_sample_names, drop = FALSE]
if (!all(colnames(filtered_counts) == sorted_sample_names)) {
  stop("Error: Column names in filtered_counts do not match sorted sample names. Possible mismatch detected.")
}
write_gct <- function(data_matrix, output_file) {
  num_genes <- nrow(data_matrix)
  num_samples <- ncol(data_matrix)
  fileConn <- file(output_file, "w")
  writeLines("#1.2", fileConn)
  writeLines(paste(num_genes, num_samples, sep="\t"), fileConn)
  gene_ids <- rownames(data_matrix)
  description <- ifelse(is.na(gene_ids) | gene_ids == "", "Unknown", "NA")  # Handle missing descriptions
  gct_data <- cbind(Name=gene_ids, Description=description, data_matrix)
  write.table(gct_data, file=fileConn, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  close(fileConn)
}
output_file <- paste0("deseq_normalized_counts_for_", args$chemical, args$conc, ".gct")
write_gct(filtered_counts, output_file)
message("GCT file successfully written: ", output_file)

sink(file = log_file, append = TRUE)
write_res <- results(dds)
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(write_res),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
gene_symbols <- as.character(gene_symbols)
gene_symbols[is.na(gene_symbols)] <- rownames(write_res)[is.na(gene_symbols)]
stopifnot(length(gene_symbols) == nrow(write_res)) # check
df <- data.frame(
  gene = gene_symbols,           # Gene symbols
  ensembl_id = rownames(write_res),  #  Ensembl IDs
  stat = write_res$padj,         # Adjusted p-values
  log2FoldChange = write_res$log2FoldChange,  # Log2 fold change
  pvalue = write_res$pvalue,     # Raw p-values
  baseMean = write_res$baseMean,  # Mean of normalized counts
  stringsAsFactors = FALSE       # ?
)
df <- df %>% filter(!is.na(stat) & !is.na(log2FoldChange) & !is.na(pvalue) & !is.na(baseMean))
df_sorted <- df[order(df$stat, na.last = TRUE), ]
write.csv(df_sorted, "fileFORbiostatsquidALLSTAT.csv", row.names = FALSE)
cat("File successfully saved: fileFORbiostatsquidALLSTAT.csv\n")

gene_symbols[is.na(gene_symbols)] <- rownames(write_res)[is.na(gene_symbols)]
df <- data.frame(
  gene = gene_symbols,           # Gene symbols
  ensembl_id = rownames(write_res),  # Ensembl IDs for reference
  stat = write_res$padj,         # Adjusted p-values
  log2FoldChange = write_res$log2FoldChange,  # Log2 fold change
  stringsAsFactors = FALSE       # ?
)
df <- df %>% filter(!is.na(stat) & !is.na(log2FoldChange))
df_sorted <- df[order(df$stat, na.last = TRUE), ]
significant_genes <- df_sorted[which(!is.na(df_sorted$stat) & df_sorted$stat < 0.05), ]
write.csv(significant_genes, paste0("sigGenes_", args$chemical, args$conc, ".csv"), row.names = FALSE)
