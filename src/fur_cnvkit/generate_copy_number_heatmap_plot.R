#!/usr/bin/env Rscript

# ----------------------------
# Load necessary libraries
# ----------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(randomcoloR)
})

# ----------------------------
# Define functions for modularity
# ----------------------------

# 1. Function to parse command-line options
parse_options <- function() {
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
    make_option(c("-c", "--chromosome_mapping"), type = "character",
                help = "Input file mapping genes to chromosomes (tab-delimited with columns: Gene and Chromosome)", metavar = "FILE"),
    make_option(c("-s", "--summary"), type = "character", default = NULL,
                help = "Output file for overall summary statistics [default: derived from OUTPUT]", metavar = "FILE"),
    make_option(c("-g", "--gene_stats"), type = "character", default = NULL,
                help = "Output file for per gene statistics [default: derived from OUTPUT]", metavar = "FILE"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Print verbose messages during processing")
  )

  opt_parser <- OptionParser(option_list = option_list)
  opts <- parse_args(opt_parser)

  # Check required arguments
  if (is.null(opts$input) || is.null(opts$output) || is.null(opts$pdf_width) ||
      is.null(opts$pdf_height) || is.null(opts$min_alterations) || is.null(opts$chromosome_mapping)) {
    print_help(opt_parser)
    stop("Missing one or more required arguments.", call. = FALSE)
  }

  # Derive default filenames if not provided
  if (is.null(opts$summary)) {
    opts$summary <- sub("\\.pdf$", "_summary.txt", opts$output)
  }
  if (is.null(opts$gene_stats)) {
    opts$gene_stats <- sub("\\.pdf$", "_gene_stats.txt", opts$output)
  }

  return(opts)
}

# 2. Function to read and filter the copy number data
read_and_filter_data <- function(input_file, min_alterations, verbose = FALSE) {
  if (verbose) cat("Reading input CSV file...\n")
  df <- read.csv(input_file, header = TRUE, row.names = 1)

  mat <- as.matrix(df)
  mat[is.na(mat)] <- 0

  # Transpose so rows are genes and columns are samples
  mat <- t(mat)

  if (verbose) cat("Filtering rows and columns with no alterations...\n")
  mat <- mat[, colSums(abs(mat)) > 0]
  mat <- mat[rowSums(abs(mat)) > 0, ]

  if (verbose) cat("Applying minimum alteration threshold...\n")
  mat <- mat[rowSums(abs(mat) > 0) >= min_alterations, ]
  mat <- mat[, colSums(abs(mat)) > 0]

  return(mat)
}

# 3. Function to prepare chromosome mapping and annotation
prepare_annotation <- function(mapping_file, gene_names, verbose = FALSE) {
  if (verbose) cat("Reading chromosome mapping file...\n")
  chromosomes <- read.table(mapping_file, header = FALSE, sep = "\t",
                            col.names = c("Gene", "Chromosome"))
  chromosomes$Gene <- as.character(chromosomes$Gene)
  chromosomes$Chromosome <- as.factor(chromosomes$Chromosome)

  # Warn if some genes are missing in the mapping file
  missing_genes <- setdiff(gene_names, chromosomes$Gene)
  if (length(missing_genes) > 0) {
    warning("The following genes are missing in the mapping file and will be assigned as 'Unknown': ",
            paste(missing_genes, collapse = ", "))
    missing_df <- data.frame(Gene = missing_genes, Chromosome = factor("Unknown"))
    chromosomes <- rbind(chromosomes, missing_df)
  }

  # Assign unique colors using randomcoloR
  unique_chromosomes <- sort(unique(chromosomes$Chromosome))
  n_chromosomes <- length(unique_chromosomes)
  pal <- distinctColorPalette(n_chromosomes)  # Generate distinct colors
  chromosome_colors <- setNames(pal, unique_chromosomes)

  # Create row annotation for the heatmap
  annotation_df <- data.frame(Chromosome = chromosomes$Chromosome[match(gene_names, chromosomes$Gene)])
  gene_annotation <- HeatmapAnnotation(
    df = annotation_df,
    col = list(Chromosome = chromosome_colors),
    which = "row"
  )
  return(gene_annotation)
}

# 4. Function to generate and save the heatmap PDF
generate_heatmap <- function(mat, gene_annotation, output_pdf, pdf_width, pdf_height, verbose = FALSE) {
  max_val <- max(abs(mat), na.rm = TRUE)
  if (verbose) cat("Generating heatmap...\n")

  heatmap_obj <- Heatmap(
    mat,
    name = "Copy Number Change",
    col = colorRamp2(c(-max_val, 0, max_val), c("blue", "white", "red")),
    left_annotation = gene_annotation,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    column_names_side = "bottom",
    heatmap_legend_param = list(
      title = "Copy Number Change",
      at = c(-max_val, 0, max_val),
      labels = c(paste0("-", max_val), "0", max_val)
    )
  )

  if (verbose) cat("Saving heatmap to PDF...\n")
  pdf(output_pdf, width = pdf_width, height = pdf_height)
  draw(heatmap_obj)
  dev.off()
}

# 5. Function to compute overall statistics
compute_overall_stats <- function(mat) {
  genes_alteration_count <- rowSums(mat != 0)
  samples_alteration_count <- colSums(mat != 0)

  # Overall stats
  overall_stats <- list(
    genes_with_alteration = sum(genes_alteration_count > 0),
    samples_with_alteration = sum(samples_alteration_count > 0),
    max_alterations_gene = rownames(mat)[which.max(genes_alteration_count)],
    max_alterations_value = max(genes_alteration_count),
    max_alterations_sample = colnames(mat)[which.max(samples_alteration_count)],
    max_sample_alterations_value = max(samples_alteration_count),
    top_genes = sort(genes_alteration_count, decreasing = TRUE)[1:12],
    top_samples = sort(samples_alteration_count, decreasing = TRUE)[1:12],
    gains_count = rowSums(mat > 0),
    losses_count = rowSums(mat < 0)
  )

  # Gains and losses
  overall_stats$total_gains <- sum(overall_stats$gains_count)
  overall_stats$total_losses <- sum(overall_stats$losses_count)
  overall_stats$avg_gains <- mean(overall_stats$gains_count)
  overall_stats$avg_losses <- mean(overall_stats$losses_count)
  overall_stats$median_gains <- median(overall_stats$gains_count)
  overall_stats$median_losses <- median(overall_stats$losses_count)
  overall_stats$gene_max_gains <- rownames(mat)[which.max(overall_stats$gains_count)]
  overall_stats$max_gains_value <- max(overall_stats$gains_count)
  overall_stats$gene_max_losses <- rownames(mat)[which.max(overall_stats$losses_count)]
  overall_stats$max_losses_value <- max(overall_stats$losses_count)
  overall_stats$top_gains_genes <- sort(overall_stats$gains_count, decreasing = TRUE)[1:12]
  overall_stats$top_losses_genes <- sort(overall_stats$losses_count, decreasing = TRUE)[1:12]

  return(overall_stats)
}

# 6. Function to write overall summary to file
write_overall_summary <- function(stats, summary_file) {
  overall_conn <- file(summary_file, "w")

  cat("----- Overall Summary -----\n", file = overall_conn)
  cat("Number of genes with a copy number alteration:", stats$genes_with_alteration, "\n", file = overall_conn)
  cat("Number of samples with a copy number alteration:", stats$samples_with_alteration, "\n", file = overall_conn)
  cat("Gene with the most alterations:", stats$max_alterations_gene, "with", stats$max_alterations_value, "alterations\n", file = overall_conn)
  cat("Sample with the most alterations:", stats$max_alterations_sample, "with", stats$max_sample_alterations_value, "alterations\n", file = overall_conn)

  cat("\n--- Top 12 Genes with the Most Alterations ---\n", file = overall_conn)
  for (i in seq_along(stats$top_genes)) {
    gene_name <- names(stats$top_genes)[i]
    cat(i, ". ", gene_name, ": ", stats$top_genes[i], " alterations\n", sep = "", file = overall_conn)
  }

  cat("\n--- Top 12 Samples with the Most Alterations ---\n", file = overall_conn)
  for (i in seq_along(stats$top_samples)) {
    sample_name <- names(stats$top_samples)[i]
    cat(i, ". ", sample_name, ": ", stats$top_samples[i], " alterations\n", sep = "", file = overall_conn)
  }

  cat("\n----- Comprehensive Gains and Losses -----\n", file = overall_conn)
  cat("Total gains across all genes:", stats$total_gains, "\n", file = overall_conn)
  cat("Total losses across all genes:", stats$total_losses, "\n", file = overall_conn)
  cat("Average gains per gene:", round(stats$avg_gains, 2), "\n", file = overall_conn)
  cat("Average losses per gene:", round(stats$avg_losses, 2), "\n", file = overall_conn)
  cat("Median gains per gene:", stats$median_gains, "\n", file = overall_conn)
  cat("Median losses per gene:", stats$median_losses, "\n", file = overall_conn)
  cat("Gene with highest gains:", stats$gene_max_gains, "with", stats$max_gains_value, "gains\n", file = overall_conn)
  cat("Gene with highest losses:", stats$gene_max_losses, "with", stats$max_losses_value, "losses\n", file = overall_conn)

  cat("\n--- Top 12 Genes with the Most Gains ---\n", file = overall_conn)
  for (i in seq_along(stats$top_gains_genes)) {
    gene_name <- names(stats$top_gains_genes)[i]
    cat(i, ". ", gene_name, ": ", stats$top_gains_genes[i], " gains\n", sep = "", file = overall_conn)
  }

  cat("\n--- Top 12 Genes with the Most Losses ---\n", file = overall_conn)
  for (i in seq_along(stats$top_losses_genes)) {
    gene_name <- names(stats$top_losses_genes)[i]
    cat(i, ". ", gene_name, ": ", stats$top_losses_genes[i], " losses\n", sep = "", file = overall_conn)
  }

  close(overall_conn)
}

# 7. Function to write per-gene statistics to file
write_gene_stats <- function(mat, gene_stats_file) {
  gene_stats <- data.frame(
    Gene = rownames(mat),
    Total_Alterations = rowSums(mat != 0),
    Gains = rowSums(mat > 0),
    Losses = rowSums(mat < 0),
    Frequency = round(rowSums(mat != 0) / ncol(mat), 3)
  )

  write.table(gene_stats, file = gene_stats_file, sep = "\t", row.names = FALSE,
              quote = FALSE, col.names = TRUE)
}

# ----------------------------
# Main Script Execution
# ----------------------------
main <- function() {
  opts <- parse_options()
  if (opts$verbose) cat("Starting processing...\n")

  # Read and filter data
  mat <- read_and_filter_data(opts$input, opts$min_alterations, verbose = opts$verbose)
  if (opts$verbose) cat("Data dimensions after filtering: ", dim(mat), "\n")

  # Prepare annotation for heatmap
  gene_annotation <- prepare_annotation(opts$chromosome_mapping, rownames(mat), verbose = opts$verbose)

  # Generate and save the heatmap
  generate_heatmap(mat, gene_annotation, opts$output, opts$pdf_width, opts$pdf_height, verbose = opts$verbose)

  # Compute overall statistics
  overall_stats <- compute_overall_stats(mat)

  # Write overall summary and gene stats to separate files
  write_overall_summary(overall_stats, opts$summary)
  write_gene_stats(mat, opts$gene_stats)

  if (opts$verbose) {
    cat("Processing complete.\n")
    cat("Overall summary written to:", opts$summary, "\n")
    cat("Gene statistics written to:", opts$gene_stats, "\n")
  }
}

# Execute main if running this script directly
if (identical(environment(), globalenv())) {
  main()
}
