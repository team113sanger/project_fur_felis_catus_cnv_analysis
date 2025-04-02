#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate an oncoprint figure showing recurrent somatic mutations and CNVs across samples."
    )
    parser.add_argument(
        "maf_file",
        type=Path,
        help="Path to the MAF file containing somatic mutation information.",
    )
    parser.add_argument(
        "cnv_file", type=Path, help="Path to the CSV file containing CNV log2 values."
    )
    parser.add_argument(
        "tumour_type", type=str, help="Tumour type to appear in the plot title."
    )
    parser.add_argument(
        "--num_recurrent_samples",
        type=int,
        default=2,
        help="Minimum number of recurrent samples for an alteration to be considered (default: 2).",
    )
    parser.add_argument(
        "--plot_height",
        type=int,
        default=12,
        help="Plot height in inches (default: 12).",
    )
    parser.add_argument(
        "--plot_width", type=int, default=20, help="Plot width in inches (default: 20)."
    )
    parser.add_argument(
        "--include_option",
        type=int,
        choices=[1, 2],
        default=2,
        help="1 for genes with both alterations; 2 for either alteration or both (default: 2).",
    )
    parser.add_argument(
        "--output_file",
        type=Path,
        default=None,
        help="Path to save the plot. If not provided, the plot will be shown interactively.",
    )
    parser.add_argument(
        "--gene_chrom_file",
        type=Path,
        default=None,
        help="Optional tab-delimited file with gene location information. "
        "If provided, the file should have either two columns (Gene, Chromosome) or five columns "
        "(Gene, Ensembl, Chromosome, Start, End). Only genes found in this file are plotted.",
    )
    parser.add_argument(
        "--gene_sort",
        choices=["alterations", "position"],
        default="alterations",
        help="Sort genes on the x-axis by number of alterations (alterations) or by genomic position (position).",
    )
    return parser.parse_args()


def read_data(maf_file: Path, cnv_file: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Reads the MAF and CNV files."""
    somatic_mutations_data = pd.read_csv(maf_file, sep="\t")
    copy_number_data = pd.read_csv(cnv_file)
    return somatic_mutations_data, copy_number_data


def build_gene_chrom_dict(
    somatic_mutations_data: pd.DataFrame, gene_chrom_file: Optional[Path]
) -> Tuple[Dict[str, str], Optional[pd.DataFrame]]:
    """
    Builds a gene-to-chromosome mapping.
    If gene_chrom_file is provided and contains location information (5 columns), it returns:
      - a dictionary mapping gene to chromosome (for plotting), and
      - a DataFrame with gene location info (used for sorting by genomic position).
    Otherwise, it falls back to using the MAF file.
    """
    gene_info_df = None
    if gene_chrom_file is not None:
        try:
            # Attempt to read a file with location info: Gene, Ensembl, Chromosome, Start, End
            gene_info_df = pd.read_csv(
                gene_chrom_file,
                sep="\t",
                header=None,
                names=["Gene", "Ensembl", "Chromosome", "Start", "End"],
            )
            gene_info_df.set_index("Gene", inplace=True)
            gene_chrom_dict = gene_info_df["Chromosome"].to_dict()
        except Exception:
            # Fallback if file has only two columns: Gene, Chromosome
            gene_chrom_map = pd.read_csv(
                gene_chrom_file, sep="\t", header=None, names=["Gene", "Chromosome"]
            )
            gene_chrom_dict = gene_chrom_map.set_index("Gene")["Chromosome"].to_dict()
    else:
        gene_chrom_dict = (
            somatic_mutations_data.groupby("Hugo_Symbol")["Chromosome"]
            .first()
            .to_dict()
        )
    return gene_chrom_dict, gene_info_df


def get_recurrent_genes(
    somatic_mutations_data: pd.DataFrame,
    copy_number_data: pd.DataFrame,
    num_recurrent_samples: int,
    include_option: int,
    gene_chrom_dict: Optional[Dict[str, str]] = None,
) -> List[str]:
    """Identifies recurrently altered genes based on the provided include_option."""
    maf_gene_counts = somatic_mutations_data["Hugo_Symbol"].value_counts()
    maf_recurrent_genes = set(
        maf_gene_counts[maf_gene_counts >= num_recurrent_samples].index
    )

    cnv_gene_counts = copy_number_data.iloc[:, 1:].apply(
        lambda x: (x != 0).sum(), axis=0
    )
    cnv_recurrent_genes = set(
        cnv_gene_counts[cnv_gene_counts >= num_recurrent_samples].index
    )

    if include_option == 1:
        recurrent_genes = maf_recurrent_genes & cnv_recurrent_genes
    elif include_option == 2:
        recurrent_genes = maf_recurrent_genes | cnv_recurrent_genes
    else:
        raise ValueError(
            "Invalid value for include_option. Use 1 for both or 2 for either."
        )

    # Keep only genes present in the CNV file and, if applicable, in the gene_chrom mapping.
    recurrent_genes = [
        gene for gene in recurrent_genes if gene in copy_number_data.columns
    ]
    if gene_chrom_dict is not None:
        recurrent_genes = [gene for gene in recurrent_genes if gene in gene_chrom_dict]
    return recurrent_genes


def prepare_cnv_heatmap_data(
    copy_number_data: pd.DataFrame, recurrent_genes: List[str]
) -> pd.DataFrame:
    """Extracts and processes CNV data for recurrent genes."""
    cna_data = copy_number_data[["Unnamed: 0"] + recurrent_genes].set_index(
        "Unnamed: 0"
    )
    return cna_data.applymap(lambda x: 1 if x > 0 else (-1 if x < 0 else 0))


def prepare_mutations_data(
    somatic_mutations_data: pd.DataFrame, recurrent_genes: List[str]
) -> pd.DataFrame:
    """Pivots the somatic mutation data for recurrent genes."""
    mutations_filtered = somatic_mutations_data[
        somatic_mutations_data["Hugo_Symbol"].isin(recurrent_genes)
    ]
    mutations_pivot = mutations_filtered.pivot_table(
        index="Tumor_Sample_Barcode",
        columns="Hugo_Symbol",
        aggfunc="size",
        fill_value=0,
    )
    return mutations_pivot.reindex(columns=recurrent_genes, fill_value=0)


def align_samples(
    cna_heatmap_data: pd.DataFrame, mutations_pivot: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    """Ensures that both CNV and mutation data are aligned to the same sample order."""
    all_samples = sorted(set(cna_heatmap_data.index).union(mutations_pivot.index))
    aligned_cna = cna_heatmap_data.reindex(all_samples).fillna(0).astype(int)
    aligned_mutations = mutations_pivot.reindex(all_samples).fillna(0).astype(int)
    return aligned_cna, aligned_mutations, all_samples


def combine_cell_value(cna_value: int, mutation_count: int) -> int:
    """
    Combine CNA value with mutation count.
    Returns:
      3  if gain (1) and mutation (>0),
     -3  if loss (-1) and mutation (>0),
      2  if neutral (0) and mutation (>0),
      or the original CNA value if no mutation.
    """
    if mutation_count > 0:
        if cna_value == 1:
            return 3
        elif cna_value == -1:
            return -3
        else:
            return 2
    return cna_value


def combine_alterations(
    aligned_cna: pd.DataFrame, aligned_mutations: pd.DataFrame
) -> pd.DataFrame:
    """Combine CNA and mutation data cell by cell."""
    combined = aligned_cna.copy()
    for sample in aligned_mutations.index:
        for gene in aligned_mutations.columns:
            mutation_count = aligned_mutations.loc[sample, gene]
            combined.loc[sample, gene] = combine_cell_value(
                combined.loc[sample, gene], mutation_count
            )
    return combined


# ----------------- Plotting Helpers ----------------- #


def draw_alteration_cell(ax, gene_index: int, sample_index: int, value: int) -> None:
    """
    Draws a single cell on the oncoprint.
    Uses a mapping: 1 -> CNA gain (lightpink), -1 -> CNA loss (lightblue),
    2 -> mutation only (star), 3 -> gain + mutation, -3 -> loss + mutation.
    """
    mapping = {
        1: ("lightpink", False),
        -1: ("lightblue", False),
        2: (None, True),
        3: ("lightpink", True),
        -3: ("lightblue", True),
    }
    rect_color, marker_flag = mapping.get(value, (None, False))
    if rect_color is not None:
        ax.add_patch(
            mpatches.Rectangle(
                (gene_index, sample_index), 1, 1, color=rect_color, zorder=1
            )
        )
    if marker_flag:
        ax.plot(
            gene_index + 0.5,
            sample_index + 0.5,
            marker="*",
            color="red",
            markersize=10,
            zorder=2,
        )


def draw_alteration_grid(ax, ordered_data: pd.DataFrame) -> None:
    """Loops over the data grid and draws each cell."""
    for sample_index, sample in enumerate(ordered_data.index):
        for gene_index, gene in enumerate(ordered_data.columns):
            value = ordered_data.loc[sample, gene]
            draw_alteration_cell(ax, gene_index, sample_index, value)


def setup_main_axis(ax, ordered_data: pd.DataFrame) -> None:
    """Configures the main axis with ticks, labels, and limits."""
    ax.set_xticks(np.arange(len(ordered_data.columns)) + 0.5)
    ax.set_xticklabels(ordered_data.columns, rotation=90)
    ax.set_yticks(np.arange(len(ordered_data.index)) + 0.5)
    ax.set_yticklabels(ordered_data.index)
    ax.set_xlim(0, len(ordered_data.columns))
    ax.set_ylim(0, len(ordered_data.index))
    ax.invert_yaxis()


def add_chromosome_twin_axis(
    ax, ordered_data: pd.DataFrame, gene_chrom_dict: Dict[str, str]
) -> Dict[str, any]:
    """
    Adds a twin x-axis to annotate chromosomes.
    Returns a dictionary with the twin axis, chrom_colors, and unique chromosomes.
    """
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    chrom_labels = [gene_chrom_dict.get(gene, "NA") for gene in ordered_data.columns]
    unique_chroms = sorted(set(chrom_labels))
    n = len(unique_chroms)
    chrom_colors = {chrom: plt.cm.tab20(i / n) for i, chrom in enumerate(unique_chroms)}
    ax2.set_xticks(np.arange(len(ordered_data.columns)) + 0.5)
    ax2.set_xticklabels([" "] * len(chrom_labels), rotation=90)
    ax2.xaxis.set_ticks_position("top")
    for label, chrom in zip(ax2.get_xticklabels(), chrom_labels):
        label.set_bbox(
            dict(
                facecolor=chrom_colors[chrom],
                alpha=1,
                edgecolor="none",
                boxstyle="square,pad=0.2",
            )
        )
    return {"ax2": ax2, "chrom_colors": chrom_colors, "unique_chroms": unique_chroms}


def add_legends(ax, chrom_colors: dict, unique_chroms: list) -> None:
    """Adds legends for CNA alterations and chromosomes."""
    gain_patch = mpatches.Patch(color="lightpink", label="CNA Gain")
    loss_patch = mpatches.Patch(color="lightblue", label="CNA Loss")
    mutation_marker = plt.Line2D(
        [0],
        [0],
        marker="*",
        color="w",
        label="Somatic Mutation",
        markerfacecolor="red",
        markersize=10,
    )
    main_legend = ax.legend(
        handles=[gain_patch, loss_patch, mutation_marker],
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0.0,
    )
    ax.add_artist(main_legend)
    chrom_patches = [
        mpatches.Patch(color=chrom_colors[ch], label=ch) for ch in unique_chroms
    ]
    ax.legend(
        handles=chrom_patches,
        title="Chromosomes",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.7),
        borderaxespad=0.0,
    )


def plot_oncoprint(
    ordered_data: pd.DataFrame,
    gene_chrom_dict: Dict[str, str],
    tumour_type: str,
    plot_height: int,
    plot_width: int,
    output_file: Optional[Path],
    x_label: str,
) -> None:
    """Creates and displays (or saves) the oncoprint plot using helper functions."""
    fig, ax = plt.subplots(figsize=(plot_width, plot_height))
    draw_alteration_grid(ax, ordered_data)
    setup_main_axis(ax, ordered_data)
    twin_info = add_chromosome_twin_axis(ax, ordered_data, gene_chrom_dict)
    plt.subplots_adjust(right=0.8)
    add_legends(ax, twin_info["chrom_colors"], twin_info["unique_chroms"])
    ax.set_title(
        f"Genetic Alterations (Somatic Mutations and Copy Number Variations) in {tumour_type}"
    )
    ax.set_xlabel(x_label)
    ax.set_ylabel("Samples")

    if output_file:
        plt.savefig(output_file, bbox_inches="tight")
    else:
        plt.show()


def order_genes_by_alterations(filtered_data: pd.DataFrame) -> pd.DataFrame:
    """Orders genes by the total number of alterations."""
    alteration_counts = (filtered_data != 0).sum(axis=0)
    ordered_genes = alteration_counts.sort_values(ascending=False).index
    return filtered_data[ordered_genes]


def filter_recurrent_genes_post_alignment(
    combined_heatmap_data: pd.DataFrame,
    aligned_mutations: pd.DataFrame,
    aligned_cna: pd.DataFrame,
    num_recurrent_samples: int,
    include_option: int,
) -> pd.DataFrame:
    """Post-process filtering of genes based on recurrent alterations across samples."""
    if include_option == 1:
        recurrent_mutation_genes = combined_heatmap_data.columns[
            (aligned_mutations > 0).sum(axis=0) >= num_recurrent_samples
        ]
        recurrent_cnv_genes = combined_heatmap_data.columns[
            (aligned_cna.abs() > 0).sum(axis=0) >= num_recurrent_samples
        ]
        recurrent_genes = list(set(recurrent_mutation_genes) & set(recurrent_cnv_genes))
    else:
        recurrent_genes = combined_heatmap_data.columns[
            (combined_heatmap_data != 0).sum(axis=0) >= num_recurrent_samples
        ]
    return combined_heatmap_data[recurrent_genes]


def generate_recurrent_gene_oncoprint(
    maf_file: Path,
    cnv_file: Path,
    tumour_type: str,
    num_recurrent_samples: int = 2,
    plot_height: int = 12,
    plot_width: int = 20,
    include_option: int = 2,
    output_file: Optional[Path] = None,
    gene_chrom_file: Optional[Path] = None,
    gene_sort: str = "alterations",
) -> None:
    # Step 1: Read data
    somatic_mutations_data, copy_number_data = read_data(maf_file, cnv_file)

    # Step 1a: Build gene-to-chromosome mapping and gene location info.
    gene_chrom_dict, gene_info_df = build_gene_chrom_dict(
        somatic_mutations_data, gene_chrom_file
    )

    # Step 2 & 3: Identify recurrent genes.
    recurrent_genes = get_recurrent_genes(
        somatic_mutations_data,
        copy_number_data,
        num_recurrent_samples,
        include_option,
        gene_chrom_dict if gene_chrom_file is not None else None,
    )

    # Step 4: Prepare CNV and mutation data.
    cna_heatmap_data = prepare_cnv_heatmap_data(copy_number_data, recurrent_genes)
    mutations_pivot = prepare_mutations_data(somatic_mutations_data, recurrent_genes)

    # Step 5: Align sample order across datasets.
    aligned_cna, aligned_mutations, _ = align_samples(cna_heatmap_data, mutations_pivot)

    # Step 6: Combine the CNA and mutation data.
    combined_heatmap_data = combine_alterations(aligned_cna, aligned_mutations)

    # Step 7: Filter genes based on recurrence post-alignment.
    filtered_data = filter_recurrent_genes_post_alignment(
        combined_heatmap_data,
        aligned_mutations,
        aligned_cna,
        num_recurrent_samples,
        include_option,
    )
    # Step 8: Order genes.
    ordered_data = order_genes_by_alterations(filtered_data)

    # If gene_sort is 'position', reorder genes based on genomic coordinates.
    if gene_sort == "position":
        if gene_info_df is None:
            raise ValueError(
                "Genomic position sorting requested but no gene_chrom_file with location info provided."
            )
        # Subset gene_info_df to the genes present in ordered_data.
        gene_info_sub = gene_info_df.loc[ordered_data.columns]
        # Sort by Chromosome and then by Start coordinate.
        sorted_genes = gene_info_sub.sort_values(
            by=["Chromosome", "Start"]
        ).index.tolist()
        ordered_data = ordered_data[sorted_genes]
        x_label = "Genes (Ordered by Genomic Position)"
    else:
        x_label = "Genes (Ordered by Number of Alterations)"

    # Step 9: Plot the oncoprint.
    plot_oncoprint(
        ordered_data,
        gene_chrom_dict,
        tumour_type,
        plot_height,
        plot_width,
        output_file,
        x_label,
    )


def main():
    args = parse_args()
    generate_recurrent_gene_oncoprint(
        maf_file=args.maf_file,
        cnv_file=args.cnv_file,
        tumour_type=args.tumour_type,
        num_recurrent_samples=args.num_recurrent_samples,
        plot_height=args.plot_height,
        plot_width=args.plot_width,
        include_option=args.include_option,
        output_file=args.output_file,
        gene_chrom_file=args.gene_chrom_file,
        gene_sort=args.gene_sort,
    )


if __name__ == "__main__":
    main()
