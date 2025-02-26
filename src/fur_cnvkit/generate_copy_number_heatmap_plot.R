#!/usr/bin/env Rscript

# Load necessary libraries
library(optparse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Input CSV file containing copy number alteration data", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output PDF file for the heatmap", metavar = "FILE"),
  make_option(c("-w", "--pdf_width"), type = "numeric",
              help = "Width of the output PDF", metavar = "WIDTH"),
  make_option(c("-t", "--pdf_height"), type = "numeric",
              help = "Height of the output PDF", metavar = "HEIGHT"),
  make_option(c("-m", "--min_alterations"), type = "integer",
              help = "Minimum number of alterations a gene must have to be included", metavar = "MIN"),
  make_option(c("-s", "--species"), type = "character",
              help = "Species: choose either 'cat' or 'dog'", metavar = "SPECIES")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Check if all required options are provided
if (is.null(opts$input) || is.null(opts$output) || is.null(opts$pdf_width) ||
    is.null(opts$pdf_height) || is.null(opts$min_alterations) || is.null(opts$species)) {
  print_help(opt_parser)
  stop("Missing one or more required arguments.", call. = FALSE)
}

# Generate summary file name by replacing .pdf with _summary.txt in the output file name
summary_file <- sub("\\.pdf$", "_summary.txt", opts$output)

# Read the CSV file
df <- read.csv(opts$input, header = TRUE, row.names = 1)

# Convert the dataframe to a matrix (ComplexHeatmap works with matrix format)
mat <- as.matrix(df)
mat[is.na(mat)] <- 0

# Transpose the matrix to flip the axes (genes become rows and samples become columns)
mat <- t(mat)

# Filter matrix to exclude rows/columns with no significant hits
mat <- mat[, colSums(abs(mat)) > 0]
mat <- mat[rowSums(abs(mat)) > 0, ]

# Filter to include only recurrently altered genes based on min_alterations
mat <- mat[rowSums(abs(mat) > 0) >= opts$min_alterations, ]

# Reapply column filter in case some columns are all zeros after row filtering
mat <- mat[, colSums(abs(mat)) > 0]

# Set the chromosome mapping file and color palette based on species
if (opts$species == "cat") {
  chromosome_file <- '/lustre/scratch125/casm/team113da/users/bf14/variant_caller_benchmarking/CNVkit/fur_cnvkit/metadata/Felis_catus_9.0.gene_chromosome_mapping.txt'
  chromosome_colors <- c(A1 = "red", A2 = "blue", A3 = "lightgreen",
                         B1 = "yellow", B2 = "orange", B3 = "purple",
                         B4 = "brown", C1 = "pink", C2 = "cyan",
                         D1 = "magenta", D2 = "grey", D3 = "darkolivegreen",
                         D4 = "limegreen", E1 = "navy", E2 = "darkcyan",
                         E3 = "maroon", F1 = "violet", F2 = "tan", X = "black")
} else if (opts$species == "dog") {
  chromosome_file <- '/lustre/scratch125/casm/team113da/users/bf14/variant_caller_benchmarking/CNVkit/fur_cnvkit/metadata/CanFam3.1.gene_chromosome_mapping.txt'
  chromosome_colors <- setNames(
    colorRampPalette(c("red", "blue", "green", "purple", "orange", "pink", "cyan", "magenta", "yellow", "brown", "navy", "grey"))(39),
    c(as.character(1:38), "X")
  )
} else {
  stop("Invalid species. Choose either 'cat' or 'dog'.")
}

# Read the chromosome mapping file
chromosomes <- read.table(chromosome_file, header = FALSE, sep = '\t',
                          col.names = c("Gene", "Chromosome"))
chromosomes$Gene <- as.character(chromosomes$Gene)
chromosomes$Chromosome <- as.factor(chromosomes$Chromosome)

# Create a row annotation for the heatmap using the chromosome mapping
gene_annotation <- HeatmapAnnotation(
  df = data.frame(Chromosome = chromosomes$Chromosome[match(rownames(mat), chromosomes$Gene)]),
  col = list(Chromosome = chromosome_colors),
  which = "row"
)

# Create the heatmap with the specified options and annotation
heatmap <- Heatmap(
  mat,
  name = "Copy Number Change",
  col = colorRamp2(c(-max(abs(mat), na.rm = TRUE), 0, max(abs(mat), na.rm = TRUE)),
                   c("blue", "white", "red")),
  left_annotation = gene_annotation,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  heatmap_legend_param = list(
    title = "Copy Number Change",
    at = c(-max(abs(mat), na.rm = TRUE), 0, max(abs(mat), na.rm = TRUE)),
    labels = c(paste0("-", max(abs(mat), na.rm = TRUE)), "0", max(abs(mat), na.rm = TRUE))
  )
)

# Save the heatmap to a PDF file
pdf(opts$output, width = opts$pdf_width, height = opts$pdf_height)
draw(heatmap)
dev.off()

# Calculate various statistics
genes_with_alteration <- sum(rowSums(mat != 0) > 0)
samples_with_alteration <- sum(colSums(mat != 0) > 0)

genes_alteration_count <- rowSums(mat != 0)
max_alterations_gene <- rownames(mat)[which.max(genes_alteration_count)]

samples_alteration_count <- colSums(mat != 0)
max_alterations_sample <- colnames(mat)[which.max(samples_alteration_count)]

# Get the top 12 genes and samples with the most alterations
top_genes <- sort(genes_alteration_count, decreasing = TRUE)[1:12]
top_samples <- sort(samples_alteration_count, decreasing = TRUE)[1:12]

# Calculate the number of gains and losses for each gene
gains_count <- rowSums(mat > 0)
losses_count <- rowSums(mat < 0)

# Identify the top 12 genes with the most gains and losses
top_gains_genes <- sort(gains_count, decreasing = TRUE)[1:12]
top_losses_genes <- sort(losses_count, decreasing = TRUE)[1:12]

# Write the summary statistics to a text file
file_conn <- file(summary_file, "w")

cat("Number of genes with a copy number alteration:", genes_with_alteration, "\n", file = file_conn)
cat("Number of samples with a copy number alteration:", samples_with_alteration, "\n", file = file_conn)
cat("Gene with the most observed copy number alterations:", max_alterations_gene, "with", max(genes_alteration_count), "alterations\n", file = file_conn)
cat("Sample with the most observed copy number alterations:", max_alterations_sample, "with", max(samples_alteration_count), "alterations\n", file = file_conn)

cat("\nTop 12 Genes with the Most Alterations:\n", file = file_conn)
for (i in seq_along(top_genes)) {
  gene_name <- names(top_genes)[i]
  alterations_count <- top_genes[i]
  cat(i, ". ", gene_name, ": ", alterations_count, " alterations\n", sep = "", file = file_conn)
}

cat("\nTop 12 Samples with the Most Alterations:\n", file = file_conn)
for (i in seq_along(top_samples)) {
  sample_name <- names(top_samples)[i]
  alterations_count <- top_samples[i]
  cat(i, ". ", sample_name, ": ", alterations_count, " alterations\n", sep = "", file = file_conn)
}

cat("\nTop 12 Genes with the Most Gains:\n", file = file_conn)
for (i in seq_along(top_gains_genes)) {
  gene_name <- names(top_gains_genes)[i]
  gains <- top_gains_genes[i]
  cat(i, ". ", gene_name, ": ", gains, " gains\n", sep = "", file = file_conn)
}

cat("\nTop 12 Genes with the Most Losses:\n", file = file_conn)
for (i in seq_along(top_losses_genes)) {
  gene_name <- names(top_losses_genes)[i]
  losses <- top_losses_genes[i]
  cat(i, ". ", gene_name, ": ", losses, " losses\n", sep = "", file = file_conn)
}

close(file_conn)
