import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

MAF_PATH = "data/tcga_laml.maf"
RANDOM_STATE = 42

DRIVER_GENES = {
    "FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2", "TET2", "RUNX1", "TP53", "NRAS", "KRAS",
    "CEBPA", "WT1", "ASXL1", "KIT", "PTPN11", "U2AF1", "SRSF2", "SF3B1", "STAG2", "PHF6",
    "EZH2", "BCOR", "GATA2", "CBL", "JAK2",
}


def strip_to_raw_variant_caller_output(maf):
    """What a variant caller actually gives you: position + alleles. No gene, no consequence."""
    raw = maf[[
        "Chromosome", "Start_Position", "End_position",
        "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
    ]].copy()
    raw.columns = ["chrom", "start", "end", "ref", "alt", "sample"]
    return raw


def build_gene_reference(train_maf):
    """
    Build a (chrom, gene) -> (min_pos, max_pos) lookup table from the training split.
    This is a simplification of a real transcript model (which uses exact exon
    boundaries) -- here it's just the observed span of mutation positions per gene
    in the training data.
    """
    ref = (
        train_maf.groupby(["Chromosome", "Hugo_Symbol"])["Start_Position"]
        .agg(["min", "max", "count"])
        .reset_index()
        .rename(columns={"min": "start", "max": "end", "count": "n_train_mutations"})
    )
    # pad each gene's span slightly, since a handful of training mutations
    # will rarely span the gene's true full length
    pad = 2000
    ref["start"] = ref["start"] - pad
    ref["end"] = ref["end"] + pad
    return ref


def annotate_gene(raw_variant, gene_ref):
    """Look up which gene (if any) a variant's position falls inside, per chromosome."""
    candidates = gene_ref[
        (gene_ref["Chromosome"] == raw_variant["chrom"]) &
        (gene_ref["start"] <= raw_variant["start"]) &
        (gene_ref["end"] >= raw_variant["start"])
    ]
    if len(candidates) == 0:
        return None, 0
    if len(candidates) == 1:
        return candidates.iloc[0]["Hugo_Symbol"], 1
    # ambiguous: more than one gene's span covers this position -- pick the
    # one whose training mutation count is highest (most likely real gene)
    best = candidates.sort_values("n_train_mutations", ascending=False).iloc[0]
    return best["Hugo_Symbol"], len(candidates)


def classify_variant_type(ref, alt):
    # MAF convention: "-" represents a zero-length allele (pure insertion/deletion),
    # not a literal one-character allele -- must be normalized before measuring length
    ref_len = 0 if ref == "-" else len(ref)
    alt_len = 0 if alt == "-" else len(alt)
    if ref_len == alt_len == 1:
        return "SNP"
    elif alt_len > ref_len:
        return "INS"
    elif alt_len < ref_len:
        return "DEL"
    else:
        return "SNP"  # multi-nucleotide substitution, same length -- treat as SNP-like


def predict_frame_consequence(ref, alt, variant_type):
    if variant_type == "SNP":
        return "n/a"
    ref_len = 0 if ref == "-" else len(ref)
    alt_len = 0 if alt == "-" else len(alt)
    length_change = abs(alt_len - ref_len)
    return "frameshift" if length_change % 3 != 0 else "in_frame"


if __name__ == "__main__":
    maf = pd.read_csv(MAF_PATH, sep="\t", low_memory=False)

    # split by SAMPLE, not by row -- so the model can't learn a gene's position
    # from one mutation in a patient and "cheat" by reusing another mutation from
    # the SAME patient at a nearby position in the test set
    samples = maf["Tumor_Sample_Barcode"].unique()
    train_samples, test_samples = train_test_split(samples, test_size=0.3, random_state=RANDOM_STATE)

    train_maf = maf[maf["Tumor_Sample_Barcode"].isin(train_samples)]
    test_maf = maf[maf["Tumor_Sample_Barcode"].isin(test_samples)]
    print(f"Training split: {len(train_maf)} mutations from {len(train_samples)} patients")
    print(f"Test split (held out): {len(test_maf)} mutations from {len(test_samples)} patients")

    gene_ref = build_gene_reference(train_maf)
    print(f"Gene reference table built: {len(gene_ref)} (chromosome, gene) entries")

    raw_test = strip_to_raw_variant_caller_output(test_maf)
    truth = test_maf[["Hugo_Symbol", "Variant_Type", "Variant_Classification"]].reset_index(drop=True)
    raw_test = raw_test.reset_index(drop=True)

    predicted_genes = []
    n_ambiguous = 0
    n_unmapped = 0
    for _, row in raw_test.iterrows():
        gene, n_candidates = annotate_gene(row, gene_ref)
        predicted_genes.append(gene)
        if n_candidates == 0:
            n_unmapped += 1
        elif n_candidates > 1:
            n_ambiguous += 1

    raw_test["predicted_gene"] = predicted_genes
    raw_test["predicted_variant_type"] = [
        classify_variant_type(r, a) for r, a in zip(raw_test["ref"], raw_test["alt"])
    ]
    raw_test["predicted_frame_consequence"] = [
        predict_frame_consequence(r, a, vt)
        for r, a, vt in zip(raw_test["ref"], raw_test["alt"], raw_test["predicted_variant_type"])
    ]
    raw_test["predicted_driver_flag"] = raw_test["predicted_gene"].isin(DRIVER_GENES)

    # --- validate against ground truth ---
    raw_test["true_gene"] = truth["Hugo_Symbol"]
    raw_test["true_variant_type"] = truth["Variant_Type"]

    covered = raw_test[raw_test["predicted_gene"].notna()]
    coverage_rate = len(covered) / len(raw_test)
    conditional_accuracy = (covered["predicted_gene"] == covered["true_gene"]).mean()
    overall_accuracy = (raw_test["predicted_gene"] == raw_test["true_gene"]).mean()
    type_correct = (raw_test["predicted_variant_type"] == raw_test["true_variant_type"]).mean()

    print(f"\n--- Validation against ground-truth MAF annotations ---")
    print(f"Reference coverage: {len(covered)} / {len(raw_test)} test variants "
          f"({coverage_rate*100:.1f}%) fell inside a known gene span")
    print(f"Accuracy WHEN a gene candidate was found: {conditional_accuracy*100:.1f}%")
    print(f"Overall accuracy (uncovered = wrong): {overall_accuracy*100:.1f}%")
    print(f"Variant type classification accuracy: {type_correct*100:.1f}%")

    driver_true = truth["Hugo_Symbol"].isin(DRIVER_GENES)
    driver_agreement = (raw_test["predicted_driver_flag"] == driver_true).mean()
    print(f"Driver-flag agreement with ground truth: {driver_agreement*100:.1f}%")

    # save a small results table
    raw_test.to_csv("annotated_output_sample.csv", index=False)
    print("\nSaved full annotated output to annotated_output_sample.csv")

    # --- plot ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    ax = axes[0]
    outcome_counts = pd.Series({
        "Correct\n(gene found)": conditional_accuracy * coverage_rate,
        "Wrong\n(gene found)": (1 - conditional_accuracy) * coverage_rate,
        "No reference\ncoverage": 1 - coverage_rate,
    })
    ax.bar(outcome_counts.index, outcome_counts.values * 100,
           color=["#3B6E8F", "#C98A3E", "#B5433E"])
    ax.set_ylabel("% of held-out test variants")
    ax.set_title("Gene Annotation Outcome\n(held-out patients)", fontweight="bold")
    for i, v in enumerate(outcome_counts.values * 100):
        ax.text(i, v + 1, f"{v:.1f}%", ha="center", fontsize=10)
    ax.set_ylim(0, 100)

    ax = axes[1]
    type_ct = pd.crosstab(raw_test["true_variant_type"], raw_test["predicted_variant_type"])
    im = ax.imshow(type_ct.values, cmap="Blues")
    ax.set_xticks(range(len(type_ct.columns))); ax.set_xticklabels(type_ct.columns)
    ax.set_yticks(range(len(type_ct.index))); ax.set_yticklabels(type_ct.index)
    ax.set_xlabel("Predicted variant type")
    ax.set_ylabel("True variant type")
    ax.set_title("Variant Type Classification", fontweight="bold")
    for i in range(len(type_ct.index)):
        for j in range(len(type_ct.columns)):
            val = type_ct.values[i, j]
            ax.text(j, i, val, ha="center", va="center",
                     color="white" if val > type_ct.values.max()/2 else "black")

    fig.tight_layout()
    fig.savefig("../../assets/img/variant_annotator.png", dpi=200, bbox_inches="tight", facecolor="white")
    print("Saved figure to ../../assets/img/variant_annotator.png")
