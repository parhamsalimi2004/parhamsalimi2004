import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

MAF_PATH = "data/tcga_laml.maf"
CLINICAL_PATH = "data/tcga_laml_annot.tsv"
GENES_OF_INTEREST = ["FLT3", "DNMT3A", "NPM1", "TP53"]

MUTANT_COLOR = "#B5433E"
WT_COLOR = "#3B6E8F"


def load_clinical(path):
    clin = pd.read_csv(path, sep="\t")
    # data quality fix: "-Inf" is used as a placeholder for missing follow-up time
    clin["days_to_last_followup"] = clin["days_to_last_followup"].replace(
        [np.inf, -np.inf], np.nan
    )
    n_before = len(clin)
    clin = clin.dropna(subset=["days_to_last_followup"])
    n_after = len(clin)
    print(f"Dropped {n_before - n_after} patients with missing follow-up time "
          f"({n_before} -> {n_after})")
    return clin


def add_mutation_status(clin, maf, gene):
    mutated_samples = set(maf.loc[maf["Hugo_Symbol"] == gene, "Tumor_Sample_Barcode"])
    clin = clin.copy()
    clin["mutant"] = clin["Tumor_Sample_Barcode"].isin(mutated_samples)
    return clin


def run_logrank(clin):
    mut = clin[clin["mutant"]]
    wt = clin[~clin["mutant"]]
    result = logrank_test(
        mut["days_to_last_followup"], wt["days_to_last_followup"],
        event_observed_A=mut["Overall_Survival_Status"],
        event_observed_B=wt["Overall_Survival_Status"],
    )
    return result.p_value, len(mut), len(wt)


def plot_panel(maf, clin, genes, out_path):
    n = len(genes)
    ncols = 2
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 4.2 * nrows))
    axes = axes.flatten()

    summary_rows = []

    for i, gene in enumerate(genes):
        ax = axes[i]
        gclin = add_mutation_status(clin, maf, gene)
        p_value, n_mut, n_wt = run_logrank(gclin)

        kmf = KaplanMeierFitter()
        kmf.fit(
            gclin.loc[gclin["mutant"], "days_to_last_followup"],
            event_observed=gclin.loc[gclin["mutant"], "Overall_Survival_Status"],
            label=f"{gene} mutant (n={n_mut})",
        )
        kmf.plot_survival_function(ax=ax, color=MUTANT_COLOR, ci_show=True)

        kmf2 = KaplanMeierFitter()
        kmf2.fit(
            gclin.loc[~gclin["mutant"], "days_to_last_followup"],
            event_observed=gclin.loc[~gclin["mutant"], "Overall_Survival_Status"],
            label=f"{gene} wild-type (n={n_wt})",
        )
        kmf2.plot_survival_function(ax=ax, color=WT_COLOR, ci_show=True)

        ax.set_title(f"{gene}  (log-rank p = {p_value:.3f})", fontweight="bold", fontsize=11)
        ax.set_xlabel("Days")
        ax.set_ylabel("Survival probability")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_ylim(0, 1.02)

        summary_rows.append({
            "gene": gene, "n_mutant": n_mut, "n_wildtype": n_wt, "p_value": round(p_value, 4)
        })

    # hide unused axes if odd number of genes
    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Overall Survival: Mutant vs. Wild-Type — TCGA-LAML", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved survival panel to {out_path}")

    return pd.DataFrame(summary_rows).sort_values("p_value")


if __name__ == "__main__":
    maf = pd.read_csv(MAF_PATH, sep="\t", low_memory=False)
    clin = load_clinical(CLINICAL_PATH)

    print(f"\nCohort with valid survival data: {len(clin)} patients")
    print(f"Events (deaths): {clin['Overall_Survival_Status'].sum()}, "
          f"Censored: {(clin['Overall_Survival_Status'] == 0).sum()}")

    summary = plot_panel(maf, clin, GENES_OF_INTEREST, "../../assets/img/survival_panel.png")
    print("\nLog-rank test summary (sorted by significance):")
    print(summary.to_string(index=False))
