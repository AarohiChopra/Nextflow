library(tidyverse)
library(fgsea)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
library(msigdbr)
library(tibble)
library(dplyr)
library(pheatmap)
library(tidyr)
library(argparse)
library(data.table)


# --- Argument Parsing ---
parser <- ArgumentParser(description = "Run FGSEA from a log2foldsorted DESeq output")
parser$add_argument("--logfile", required = TRUE, help = "Path to log2foldsorted_file_*.csv from DeseqAutomate")
parser$add_argument("--chemical", required = TRUE, help = "Chemical/treatment name (e.g., BPA)")
parser$add_argument("--conc",     required = TRUE, help = "Concentration value (e.g., 50)")
parser$add_argument("--control",  required = TRUE, help = "Control condition label (e.g., DMSO)")
parser$add_argument("--outdir",   default = ".",  help = "Output directory")

args <- parser$parse_args()

safe_label <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
args$conc     <- gsub("\\[|\\]", "", args$conc)  # strip [ ]
args$chemical <- safe_label(args$chemical)
args$conc     <- safe_label(args$conc)
args$control  <- safe_label(args$control)

# Logging
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
log_file   <- file.path(args$outdir, sprintf("fgsea_log_%s_%s_%s.txt", args$chemical, args$conc, args$control))
sink(log_file, split = TRUE)

# Make sure we always close sinks and graphics devices even on error
.on.exit <- function() {
  try(sink(NULL), silent = TRUE)
  while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
}
# ensure on-exit runs
reg.finalizer(environment(), function(...) .on.exit(), onexit = TRUE)

cat("Arguments passed in the run command:\n"); print(args)
              
plots_pdf <- file.path(args$outdir, sprintf("fgsea_plots_%s_%s_%s.pdf", args$chemical, args$conc, args$control))
pdf(plots_pdf, width = 12, height = 8, onefile = TRUE)

# ---- Inputs (keep your names) ----
control_files <- list()
key <- paste0(args$chemical, "_", args$conc)  # e.g., "BPA_50"
control_files[[key]] <- args$logfile
print(control_files)

# this method genes the gene pathways and genes associated to the subcomponent and then converts it into a dictionary of lists like structure
get_pathway_genes <- function(subcat){
  data <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = subcat)
  subgroups <- split(data$gene_symbol, data$gs_name)
  return(subgroups)
}
# this methods is just converting the ensemble IDs from featureCounts to gene symbols and only keeping the stat and gene name column
get_stats <- function(x){
  data <- read.csv(x)
  data$gene <- gsub("\\..*", "", data$gene)
  data$gene_symbol <- mapIds(org.Hs.eg.db,
                             keys = data$gene,
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first")
  data <- na.omit(data)
  data <- data[!duplicated(data$gene_symbol), ]
  rel_cols <- data[, c("gene_symbol", "stat")]
  return(deframe(rel_cols))
}
subcategories <- c("GO:BP", "GO:CC", "GO:MF")

my_stats <- get_stats(control_files[[1]])
pathway_genes <- get_pathway_genes("GO:BP")
cat("******************************************************\n")
cat("Genes in stats:", length(my_stats), "\n")
cat("Unique genes in pathways:", length(unique(unlist(pathway_genes))), "\n") # confirming the whole package is installed
cat("Common genes:", length(intersect(names(my_stats), unique(unlist(pathway_genes)))), "\n")

sig_pathway_counts <- list() # for counts of pathways (this is a list of dataframes)

sig_pathways_by_component <- list( # for common pathways
  "GO:BP" = list(),
  "GO:CC" = list(),
  "GO:MF" = list()
)

heatmap_data_list <- list() # for heatmaps

for (i in names(control_files)) {
  file <- control_files[[i]]
  my_stats <- get_stats(file)
  #print(my_stats)
  my_stats
  my_stats <- my_stats[!is.na(names(my_stats)) & is.finite(my_stats)] # remove NAs and infinite values

  score_type <- if (all(my_stats > 0)) "pos" else if (all(my_stats < 0)) "neg" else "std"

  for (j in subcategories) {
    pathway_genes <- get_pathway_genes(j)
    pathway_genes <- pathway_genes[lengths(pathway_genes) > 0]
    if (length(pathway_genes) == 0) {
    warning(paste("Genes missing for:", j))
    next}
    #my_stats <- my_stats[names(my_stats) %in% unlist(pathway_genes)]

    # if (length(my_stats) == 0) {
    #   warning(paste("Something is wrong in:", j))
    #   next
    # }
    res <- fgsea(
      pathways = pathway_genes,
      stats = my_stats,
      eps = 0.0,
      scoreType = score_type,
      minSize = 10,
      maxSize = 500
    )

    sig_count <- sum(res$padj < 0.05, na.rm = TRUE)
    sig_pathway_counts[[length(sig_pathway_counts) + 1]] <- data.frame(
    concentration = i,
    category = j,
    sig_count = sig_count)

    sig_pathways <- res[padj < 0.05, pathway]
    sig_pathways_by_component[[j]][[i]] <- sig_pathways

    #top 10 up and down pathways
    res_clean <- res[!is.na(padj)]
    topPathwaysUp <- res_clean[NES > 0][order(padj)][1:10, pathway]
    topPathwaysDown <- res_clean[NES < 0][order(padj)][1:10, pathway]
    topPathways <- unique(c(topPathwaysUp, topPathwaysDown))

    top_res <- res[pathway %in% topPathways, .(pathway, NES)]
    top_res$concentration <- i
    top_res$subcomponent <- j

    heatmap_data_list[[paste(i, j, sep = "_")]] <- top_res


    plt <- plotGseaTable(
    pathway_genes[topPathways],
    my_stats,
    res,
    gseaParam = 0.5,
    colwidths = c(10, 3, 0.8, 1.2, 1.2),
    pathwayLabelStyle = list(size = 8, color = "black"),
    valueStyle = list(size = 8)
    )
    print(plt)
  }
}

#CREATING STACKED GRAPHS  FOR DIFFERENTIALLY EXPRESSED PATHWAYS:

sig_df <- bind_rows(sig_pathway_counts)
sig_df$category <- factor(sig_df$category,
                          levels = c("GO:BP", "GO:CC", "GO:MF"),
                          labels = c("Biological Process", "Cellular Component", "Molecular Function"))

keys <- unique(sig_df$concentration)
chem_prefix <- paste0("^", args$chemical, "_?")
conc_only   <- sub(chem_prefix, "", keys)             # e.g., "50"
pretty_lab  <- paste0(conc_only, " µM")               # e.g., "50 µM"
recode_map  <- stats::setNames(pretty_lab, keys)

sig_df$concentration <- recode(sig_df$concentration, !!!recode_map)
conc_num <- as.numeric(sub(" µM$", "", sig_df$concentration))
# get unique labels in numeric order
ord <- order(match(sig_df$concentration, unique(sig_df$concentration))[!duplicated(sig_df$concentration)], na.last = NA)
unique_labels <- unique(sig_df$concentration)[order(as.numeric(sub(" µM$", "", unique(sig_df$concentration))))]

sig_df$concentration <- factor(sig_df$concentration, levels = unique_labels)

g <- ggplot(sig_df, aes(x = concentration, y = sig_count, fill = category)) +
  geom_bar(stat = "identity") +
  labs(
    x = "BPA Concentration (µM)",
    y = "Number of Significantly Enriched Pathways (padj < 0.05)",
    fill = "GO Category",
    title = "Differentially Enriched GO Pathways Across BPA Doses"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(g)
#Finding commmon pathways among all concentrations:

# Safe multi-intersect: returns character(0) if <2 lists or any is NULL
common_across <- function(pathway_list, keys) {
  keys <- intersect(names(pathway_list), keys)
  if (length(keys) < 2) return(character(0))
  keep <- pathway_list[keys]
  # Drop NULL entries
  keep <- keep[!vapply(keep, is.null, logical(1))]
  if (length(keep) < 2) return(character(0))
  Reduce(intersect, keep)
}

# Available keys in your run (use any category; they share the same keys)
available_keys <- unique(unlist(lapply(sig_pathways_by_component, names)))

set1_concs <- c("5","50","10")
set1_keys  <- paste0(args$chemical, "_", set1_concs)

if (length(intersect(available_keys, set1_keys)) >= 2) {
  common_pathways <- lapply(sig_pathways_by_component, function(pl) {
    common_across(pl, set1_keys)
  })
  cat("******************************************************\n")
  cat("\nCommon (GO:BP) for 5/50/10 µM:\n");  print(common_pathways[["GO:BP"]])
  cat("\nCommon (GO:CC) for 5/50/10 µM:\n");  print(common_pathways[["GO:CC"]])
  cat("\nCommon (GO:MF) for 5/50/10 µM:\n");  print(common_pathways[["GO:MF"]])
} else {
  message("Skipping common-pathways (5/50/10 µM): fewer than 2 of these doses available in this run.")
}

# Common pathways among low doses 0.0005, 0.001, 0.01 µM (only if present) which are also closer to the estradiol conc we looked at ----
set2_concs <- c("0.0005","0.001","0.01")
set2_keys  <- paste0(args$chemical, "_", set2_concs)

if (length(intersect(available_keys, set2_keys)) >= 2) {
  common_pathways <- lapply(sig_pathways_by_component, function(pl) {
    common_across(pl, set2_keys)
  })
  cat("******************************************************\n")
  cat("\nCommon (GO:BP) for 0.0005/0.001/0.01 µM:\n"); print(common_pathways[["GO:BP"]])
  cat("\nCommon (GO:CC) for 0.0005/0.001/0.01 µM:\n"); print(common_pathways[["GO:CC"]])
  cat("\nCommon (GO:MF) for 0.0005/0.001/0.01 µM:\n"); print(common_pathways[["GO:MF"]])
} else {
  message("Skipping common-pathways (0.0005/0.001/0.01 µM): fewer than 2 of these doses available in this run.")
}

# Sum significant counts across categories per concentration
sig_df <- dplyr::bind_rows(sig_pathway_counts)
totals <- sig_df |>
  dplyr::group_by(concentration) |>
  dplyr::summarise(total = sum(sig_count, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(total))

top3_keys <- head(totals$concentration, 3)

if (length(top3_keys) >= 2) {
  common_pathways <- lapply(sig_pathways_by_component, function(pl) {
    common_across(pl, top3_keys)
  })
  cat("******************************************************\n")
  cat("\nCommon pathways across top-3 concentrations by total sig pathways:\n")
  print(common_pathways)
} else {
  message("Skipping auto top-3 common-pathways: <2 concentrations available.")
}

#COMPARING THE DATA WITH THE 50 GENE BIOMARKERS FOR ESTROGEN RESPONSE

file <- read.csv("/scratch/home/achopra/BPA_Alt_Human/BPA/Fgsea/biomarkers.csv")
biomarker_df <- file[!is.na(file$Fold.change.in.the.46.gene.ER.Biomarker), ]
biomarker_genes <- biomarker_df$Gene
reference_fc <- biomarker_df$Fold.change.in.the.46.gene.ER.Biomarker
names(reference_fc) <- biomarker_genes
biomarker_gene_set <- list("Estrogen_Biomarkers" = biomarker_genes)
lfc_matrix <- data.frame(Gene = biomarker_genes)
for (i in names(control_files)) {
  file <- control_files[[i]]
  stats <- get_stats(file)
  stats <- sort(stats, decreasing = TRUE)
  stats <- stats[!is.na(names(stats)) & is.finite(stats)]
  estrogen_res <- fgsea(pathways = biomarker_gene_set, stats = stats, scoreType = "std", minSize = 1, maxSize = 500)
  print(estrogen_res)
  topPathwaysUp <- estrogen_res[ES > 0][head(order(pval), n = 10), pathway]
  topPathwaysDown <- estrogen_res[ES < 0][head(order(pval), n = 10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  print(plotGseaTable(
    biomarker_gene_set[topPathways],
    stats,
    estrogen_res,
    gseaParam = 0.5,
    colwidths = c(10, 3, 0.8, 1.2, 1.2),
    pathwayLabelStyle = list(size = 8, color = "black"),
    valueStyle = list(size = 8)
  ))
  common_genes <- intersect(biomarker_genes, names(stats))
  estrogen_gene_stats <- data.frame(
    Gene                = common_genes,
    Log2FC              = stats[common_genes],
    Reference           = reference_fc[common_genes],
    Regulation          = ifelse(stats[common_genes] > 0, "Upregulated", "Downregulated"),
    Reference_Regulation = ifelse(reference_fc[common_genes] > 0, "Upregulated", "Downregulated")
  )
  print(estrogen_gene_stats)
  match_up <- sum(estrogen_gene_stats$Reference_Regulation == "Upregulated" & estrogen_gene_stats$Regulation == "Upregulated")
  match_down <- sum(estrogen_gene_stats$Reference_Regulation == "Downregulated" & estrogen_gene_stats$Regulation == "Downregulated")
  total_up <- sum(estrogen_gene_stats$Reference_Regulation == "Upregulated")
  total_down <- sum(estrogen_gene_stats$Reference_Regulation == "Downregulated")
  cat("******************************************************\n")
  cat("Upregulated match:", match_up, "out of", total_up, "\n")
  cat("Downregulated match:", match_down, "out of", total_down, "\n")
  chem_prefix <- paste0("^", args$chemical, "_?")
  conc_only   <- sub(chem_prefix, "", i)          # "50"
  pretty_col  <- paste0(conc_only, " µM")        # "50 µM"

  temp_df <- data.frame(Gene = common_genes, value = stats[common_genes], stringsAsFactors = FALSE)
  colnames(temp_df)[2] <- pretty_col

  lfc_matrix <- dplyr::full_join(lfc_matrix, temp_df, by = "Gene")
}

## 3) Finalize matrix (Reference first, doses after, numerically ordered)
lfc_matrix$Reference <- reference_fc[lfc_matrix$Gene]
rownames(lfc_matrix) <- lfc_matrix$Gene
lfc_matrix <- lfc_matrix[, setdiff(colnames(lfc_matrix), "Gene"), drop = FALSE]

# put "Reference" first, then order dose columns by numeric µM
dose_cols <- setdiff(colnames(lfc_matrix), "Reference")
# keep only those that end with " µM" (in case single run has exactly one dose)
dose_num  <- suppressWarnings(as.numeric(sub(" \\u00B5M$| µM$", "", dose_cols)))
dose_cols_sorted <- dose_cols[order(dose_num)]

lfc_matrix <- lfc_matrix[, c("Reference", dose_cols_sorted), drop = FALSE]

# Clean any inf/nan
lfc_matrix_mat <- as.matrix(lfc_matrix)
lfc_matrix_mat[!is.finite(lfc_matrix_mat)] <- 0

# Transpose for heatmap: rows = conditions, cols = genes
lfc_matrix_t <- t(lfc_matrix_mat)

## 4) Heatmap to outdir, named with chemical (and works for single dose)
if (nrow(lfc_matrix_t) > 0 && ncol(lfc_matrix_t) > 0) {
  # Title page
  grid::grid.newpage()
  grid::grid.text(
    sprintf("Estrogen Biomarker Response (%s, %s µM)", args$chemical, args$conc),
    gp = grid::gpar(fontsize = 10, fontface = "bold")
  )

  # Heatmap page(s)
  heat <- pheatmap(
    mat = lfc_matrix_t,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    display_numbers = round(lfc_matrix_t, 2),
    fontsize_col = 8,
    fontsize_row = 10,
    angle_col = 90,
    color = colorRampPalette(c("blue", "white", "red"))(200),
    na_col = "grey90",
    gaps_row = 1,
    border_color = NA,
    silent = TRUE      # return gtable without auto-plotting
  )
  grid::grid.newpage()
  grid::grid.draw(heat$gtable)
} else {
  message("Biomarker heatmap skipped: no overlapping genes with stats.")
}
dev.off()
sink(NULL); sink(NULL, type = "message")

## 5) Export the matrix (CSV) to outdir with clear name
out_csv <- file.path(args$outdir, sprintf("EstroBiomarker_Mat_%s_%s_.csv", args$chemical, args$conc))
write.csv(lfc_matrix_t, file = out_csv, row.names = TRUE)
