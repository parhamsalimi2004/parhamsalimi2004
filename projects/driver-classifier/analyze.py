import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix, classification_report
from sklearn.preprocessing import StandardScaler

MAF_PATH = "data/tcga_laml.maf"
RANDOM_STATE = 42

DRIVER_GENES = {
    "FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2", "TET2", "RUNX1", "TP53", "NRAS", "KRAS",
    "CEBPA", "WT1", "ASXL1", "KIT", "PTPN11", "U2AF1", "SRSF2", "SF3B1", "STAG2", "PHF6",
    "EZH2", "BCOR", "GATA2", "CBL", "JAK2",
}


def engineer_features(maf):
    df = maf.copy()
    df["label"] = df["Hugo_Symbol"].isin(DRIVER_GENES).astype(int)

    # hotspot recurrence: how many times does this exact protein change appear in the cohort?
    # (minus 1 so a mutation doesn't count itself)
    change_counts = df.groupby("Protein_Change")["Protein_Change"].transform("count")
    df["hotspot_count"] = change_counts - 1
    df["is_hotspot"] = (df["hotspot_count"] > 0).astype(int)

    df["vaf"] = df["i_TumorVAF_WU"]
    df["is_silent"] = (df["Variant_Classification"] == "Silent").astype(int)
    df["is_truncating"] = df["Variant_Classification"].isin(
        ["Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "Splice_Site"]
    ).astype(int)
    df["is_missense"] = (df["Variant_Classification"] == "Missense_Mutation").astype(int)
    df["is_indel"] = df["Variant_Type"].isin(["INS", "DEL"]).astype(int)

    feature_cols = ["vaf", "hotspot_count", "is_hotspot", "is_silent",
                     "is_truncating", "is_missense", "is_indel"]
    df = df.dropna(subset=feature_cols)
    return df, feature_cols


def train_and_evaluate(df, feature_cols):
    X = df[feature_cols]
    y = df["label"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, random_state=RANDOM_STATE, stratify=y
    )

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    models = {
        "Logistic Regression": LogisticRegression(class_weight="balanced", random_state=RANDOM_STATE),
        "Random Forest": RandomForestClassifier(
            n_estimators=300, max_depth=5, class_weight="balanced", random_state=RANDOM_STATE
        ),
    }

    results = {}
    for name, model in models.items():
        if name == "Logistic Regression":
            model.fit(X_train_scaled, y_train)
            proba = model.predict_proba(X_test_scaled)[:, 1]
        else:
            model.fit(X_train, y_train)
            proba = model.predict_proba(X_test)[:, 1]

        auc = roc_auc_score(y_test, proba)
        preds = (proba >= 0.5).astype(int)
        results[name] = {
            "model": model, "proba": proba, "auc": auc,
            "y_test": y_test, "preds": preds,
        }
        print(f"\n=== {name} ===")
        print(f"ROC-AUC: {auc:.3f}")
        print(classification_report(y_test, preds, target_names=["non-driver", "driver-gene"]))

    return results, X_train.columns.tolist()


def plot_results(results, feature_names, rf_model, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    # ROC curves
    ax = axes[0]
    for name, r in results.items():
        fpr, tpr, _ = roc_curve(r["y_test"], r["proba"])
        ax.plot(fpr, tpr, label=f"{name} (AUC={r['auc']:.3f})")
    ax.plot([0, 1], [0, 1], "--", color="gray", linewidth=1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("ROC Curves", fontweight="bold")
    ax.legend(fontsize=8)

    # confusion matrix for random forest
    ax = axes[1]
    rf_result = results["Random Forest"]
    cm = confusion_matrix(rf_result["y_test"], rf_result["preds"])
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["non-driver", "driver-gene"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["non-driver", "driver-gene"])
    ax.set_xlabel("Predicted"); ax.set_ylabel("Actual")
    ax.set_title("Random Forest — Confusion Matrix", fontweight="bold")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, cm[i, j], ha="center", va="center",
                     color="white" if cm[i, j] > cm.max() / 2 else "black", fontsize=13)

    # feature importance
    ax = axes[2]
    importances = pd.Series(rf_model.feature_importances_, index=feature_names).sort_values()
    ax.barh(importances.index, importances.values, color="#3B6E8F")
    ax.set_title("Random Forest — Feature Importance", fontweight="bold")
    ax.set_xlabel("Importance")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"\nSaved figure to {out_path}")


if __name__ == "__main__":
    maf = pd.read_csv(MAF_PATH, sep="\t", low_memory=False)
    df, feature_cols = engineer_features(maf)
    print(f"Usable mutations after dropping missing VAF: {len(df)}")
    print(f"Class balance -> driver-gene: {df['label'].sum()}, non-driver: {(df['label']==0).sum()}")

    results, feature_names = train_and_evaluate(df, feature_cols)
    plot_results(results, feature_names, results["Random Forest"]["model"],
                 "../../assets/img/driver_classifier.png")
