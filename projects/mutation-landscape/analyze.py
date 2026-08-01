import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

DATA_PATH = "data/tcga_laml.maf"
TOP_N_GENES = 20

# ---- Color convention: standard MAF Variant_Classification palette ----
# These are the same category colors used across the site (see style.css)
MUTATION_COLORS = {
    "Missense_Mutation": "#3B6E8F",
    "Nonsense_Mutation": "#B5433E",
    "Frame_Shift_Del": "#C98A3E",
    "Frame_Shift_Ins": "#D9A857",
    "Splice_Site": "#6C5B8C",
    "In_Frame_Del": "#5C8A66",
    "In_Frame_Ins": "#7FA88A",
    "Nonstop_Mutation": "#8C4B4B",
    "Translation_Start_Site": "#4B6C8C",
    "Silent": "#B8B4A8",
    "RNA": "#B8B4A8",
    "Multi_Hit": "#1B211F",
}
DEFAULT_COLOR = "#B8B4A8"


def load_data(path):
    maf = pd.read_csv(path, sep="\t", low_memory=False)
    return maf


def summary_stats(maf):
    n_samples = maf["Tumor_Sample_Barcode"].nunique()
    n_genes = maf["Hugo_Symbol"].nunique()
    n_mutations = len(maf)
    avg_mut_per_sample = n_mutations / n_samples
    return {
        "n_samples": n_samples,
        "n_genes_mutated": n_genes,
        "n_mutations": n_mutations,
        "avg_mutations_per_sample": round(avg_mut_per_sample, 2),
    }


def mutation_type_breakdown(maf):
    return maf["Variant_Classification"].value_counts()


def build_oncoplot_matrix(maf, top_n=TOP_N_GENES):
    # Genes ranked by number of samples they're mutated in (not raw mutation count)
    gene_sample_counts = (
        maf.groupby("Hugo_Symbol")["Tumor_Sample_Barcode"].nunique().sort_values(ascending=False)
    )
    top_genes = gene_sample_counts.head(top_n).index.tolist()

    samples = sorted(maf["Tumor_Sample_Barcode"].unique())

    # For each gene/sample cell, keep the variant classification (or "Multi_Hit" if >1 type)
    sub = maf[maf["Hugo_Symbol"].isin(top_genes)]
    matrix = pd.DataFrame(index=top_genes, columns=samples, dtype=object)
    matrix[:] = None

    grouped = sub.groupby(["Hugo_Symbol", "Tumor_Sample_Barcode"])["Variant_Classification"].apply(
        lambda x: x.iloc[0] if x.nunique() == 1 else "Multi_Hit"
    )
    for (gene, sample), vc in grouped.items():
        matrix.loc[gene, sample] = vc

    # Reorder samples: put samples with more mutations (among top genes) first, for a cleaner look
    mutated_counts = matrix.notna().sum(axis=0)
    ordered_samples = mutated_counts.sort_values(ascending=False).index.tolist()
    matrix = matrix[ordered_samples]

    # Reorder genes by frequency (already sorted by top_genes, but re-affirm)
    matrix = matrix.loc[top_genes]

    return matrix, gene_sample_counts.head(top_n)


def plot_oncoplot(matrix, gene_counts, n_total_samples, out_path):
    genes = matrix.index.tolist()
    samples = matrix.columns.tolist()
    n_genes = len(genes)
    n_samples = len(samples)

    fig = plt.figure(figsize=(14, 0.42 * n_genes + 2.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[6, 1], wspace=0.03)
    ax = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1], sharey=ax)

    # background grid (unmutated cells)
    ax.set_xlim(0, n_samples)
    ax.set_ylim(0, n_genes)
    ax.set_facecolor("#EDEFEA")

    for i in range(n_samples):
        for j, gene in enumerate(genes):
            val = matrix.loc[gene, samples[i]]
            if pd.notna(val):
                color = MUTATION_COLORS.get(val, DEFAULT_COLOR)
                rect = mpatches.Rectangle(
                    (i, n_genes - j - 1), 1, 1, facecolor=color, edgecolor="#EDEFEA", linewidth=0.3
                )
                ax.add_patch(rect)

    ax.set_yticks([n_genes - j - 0.5 for j in range(n_genes)])
    ax.set_yticklabels(genes, fontfamily="monospace", fontsize=10)
    ax.set_xticks([])
    ax.set_xlabel(f"{n_samples} tumor samples", fontsize=10, color="#1B211F")
    for spine in ax.spines.values():
        spine.set_visible(False)

    # side bar: % of samples mutated per gene
    pct = (gene_counts / n_total_samples * 100).values
    ax_bar.barh(
        [n_genes - j - 0.5 for j in range(n_genes)],
        pct,
        height=0.6,
        color="#1B211F",
    )
    ax_bar.set_xlabel("% samples", fontsize=9, color="#1B211F")
    ax_bar.set_yticks([])
    for spine in ax_bar.spines.values():
        spine.set_visible(False)
    ax_bar.tick_params(axis="x", labelsize=8)

    # legend
    legend_types = [t for t in MUTATION_COLORS if t != "Multi_Hit"]
    present_types = [t for t in legend_types if t in matrix.values]
    handles = [
        mpatches.Patch(color=MUTATION_COLORS[t], label=t.replace("_", " "))
        for t in present_types
    ]
    handles.append(mpatches.Patch(color=MUTATION_COLORS["Multi_Hit"], label="Multiple types"))
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
        fontsize=9,
    )

    fig.suptitle(
        f"TCGA-LAML Mutation Landscape — Top {n_genes} Genes",
        fontsize=14,
        fontweight="bold",
        color="#1B211F",
        y=1.02,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved oncoplot to {out_path}")


def plot_mutation_type_bar(breakdown, out_path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    colors = [MUTATION_COLORS.get(t, DEFAULT_COLOR) for t in breakdown.index]
    ax.barh(breakdown.index[::-1].str.replace("_", " "), breakdown.values[::-1], color=colors[::-1])
    ax.set_xlabel("Number of mutations")
    ax.set_title("Mutation Type Breakdown — TCGA-LAML", fontweight="bold")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved mutation type chart to {out_path}")


if __name__ == "__main__":
    maf = load_data(DATA_PATH)
    stats = summary_stats(maf)
    print("Summary statistics:", stats)

    breakdown = mutation_type_breakdown(maf)
    print("\nMutation type breakdown:\n", breakdown)

    matrix, gene_counts = build_oncoplot_matrix(maf, TOP_N_GENES)
    plot_oncoplot(matrix, gene_counts, stats["n_samples"], "../../assets/img/oncoplot.png")
    plot_mutation_type_bar(breakdown, "../../assets/img/mutation_types.png")

    print("\nTop mutated genes (by % of samples):")
    print((gene_counts / stats["n_samples"] * 100).round(1).head(10))
