#!/usr/bin/env python3

"""
01_prepare_decision_table.py

Prepare and evaluate a decision table for motif functionality model training.

Run as a subcommand tool. Each subcommand does one thing:

    aggregate       Collapse cell-line-specific annotation columns into a single
                    feature per annotation type (numeric: max; categorical: first).

    align           Align an aggregated feature table to a reference schema/template,
                    filling any missing columns with zeros.

    balance         Subsample positives and negatives to equal counts, then shuffle.

    clean           Drop rows with obvious label conflicts (e.g. a labeled-positive
                    row that has no open chromatin or footprint support).

    evaluate        Train a logistic regression model and report accuracy, precision,
                    recall, PR-AUC, and ROC-AUC on a held-out test split.

    compare         Plot PR and ROC curves for two prediction CSV files side by side.

    intersect-mpra  Intersect a motif BED file with MPRA region BED files using
                    pybedtools (returns motifs that fall inside MPRA-tested regions).

Example:
    python scripts/01_prepare_decision_table.py evaluate \\
        --table data/processed/decision_table.tsv \\
        --outdir results/tables/decision_table_evaluation
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    auc,
    confusion_matrix,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import train_test_split


def simplify_label(value):
    value = str(value)

    if value == "0":
        return "NO"
    if "Tss" in value or "Tx" in value or "BivFlnk" in value or "PLS" in value:
        return "TSS"
    if "Enh" in value or "ELS" in value:
        return "Enh"
    if "Repr" in value:
        return "Repr"
    if "Quies" in value or "Rpts" in value or "Het" in value:
        return "Quies"

    return value


def aggregate_features(args):
    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    prefixes = tuple(f"{x}___" for x in args.cells)

    keep_cols = [c for c in df.columns if c.startswith(prefixes)]
    if not keep_cols:
        raise ValueError(f"No columns found for cells: {', '.join(args.cells)}")

    out = df[keep_cols].copy()
    out.columns = [c.split("___")[-1].lower() for c in out.columns]

    collapsed = pd.DataFrame(index=out.index)

    for feature in sorted(set(out.columns)):
        block = out.loc[:, out.columns == feature]
        numeric = block.apply(pd.to_numeric, errors="coerce")

        if numeric.notna().any().any():
            collapsed[feature] = numeric.max(axis=1).fillna(0)
        else:
            collapsed[feature] = block.iloc[:, 0]

    for col in ["dnase__seq", "fantom", "footprints"]:
        if col in collapsed.columns:
            collapsed[col] = (
                pd.to_numeric(collapsed[col], errors="coerce")
                .fillna(0)
                .gt(0)
                .astype(int)
            )

    if "numothertfbinding" in collapsed.columns:
        collapsed["numothertfbinding"] = (
            pd.to_numeric(collapsed["numothertfbinding"], errors="coerce")
            .fillna(0)
            .clip(upper=args.max_tf_binding)
        )

    if "tfexpr" in collapsed.columns:
        x = pd.to_numeric(collapsed["tfexpr"], errors="coerce").replace(0, np.nan)
        collapsed["tfexpr"] = np.log10(x).replace([np.inf, -np.inf], 0).fillna(0)

    for state_col in ["chromhmm", "ccre"]:
        if state_col in collapsed.columns:
            simplified = collapsed[state_col].apply(simplify_label)
            dummies = pd.get_dummies(simplified)
            collapsed = pd.concat([collapsed, dummies], axis=1)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    collapsed.to_csv(args.output, sep="\t", index=False)
    print(f"Saved aggregated feature table: {args.output}")


def align_to_schema(args):
    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    template = pd.read_csv(args.template, sep="\t", nrows=1)
    schema_cols = template.columns.tolist()

    drop_cols = [
        "ccre",
        "chromhmm",
        "contactingdomain",
        "loopdomain",
        "othertfbinding",
        "replidomain",
        "tfbinding",
        "NO",
    ]

    df = df.drop(columns=[c for c in drop_cols if c in df.columns], errors="ignore")

    if args.label is not None:
        df[args.label_column] = int(args.label)

    for col in schema_cols:
        if col not in df.columns:
            df[col] = 0

    df = df[schema_cols]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Saved schema-aligned table: {args.output}")


def balance_table(args):
    pos = pd.read_csv(args.positives, sep="\t", low_memory=False)
    neg = pd.read_csv(args.negatives, sep="\t", low_memory=False)

    if args.label_column not in pos.columns or args.label_column not in neg.columns:
        raise ValueError(f"Both tables must contain '{args.label_column}'.")

    n = min(len(pos), len(neg))

    pos = pos.sample(n=n, random_state=args.seed)
    neg = neg.sample(n=n, random_state=args.seed)

    out = (
        pd.concat([pos, neg], axis=0)
        .sample(frac=1, random_state=args.seed)
        .reset_index(drop=True)
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    print(f"Saved balanced table: {args.output}")
    print(f"Rows: {len(out):,}")
    print(f"Positives: {(out[args.label_column] == 1).sum():,}")
    print(f"Negatives: {(out[args.label_column] == 0).sum():,}")


def clean_table(args):
    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    before = len(df)

    label = args.label_column

    required_cols = [label, "Repr", "Quies", "dnase__seq", "footprints", "TSS", "Enh"]
    missing = [c for c in required_cols if c not in df.columns]

    if missing:
        raise ValueError(
            "Cannot run clean step because these columns are missing: "
            + ", ".join(missing)
        )

    bad_positive = (
        (df[label] == 1)
        & (df["Repr"] == 1)
        & (df["Quies"] == 1)
        & (df["dnase__seq"] == 0)
        & (df["footprints"] == 0)
    )

    bad_negative = (
        (df[label] == 0)
        & (
            (df["TSS"] == 1)
            | ((df["Enh"] == 1) & (df["dnase__seq"] > 0))
        )
    )

    cleaned = df.loc[~(bad_positive | bad_negative)].copy()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cleaned.to_csv(args.output, sep="\t", index=False)

    print(f"Saved cleaned table: {args.output}")
    print(f"Rows before: {before:,}")
    print(f"Rows after: {len(cleaned):,}")
    print(f"Removed: {before - len(cleaned):,}")


def evaluate_table(args):
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.table, sep="\t", low_memory=False)

    if "Unnamed: 0" in df.columns:
        df = df.drop(columns=["Unnamed: 0"])

    if args.label_column not in df.columns:
        raise ValueError(f"Missing label column: {args.label_column}")

    version = args.name or args.table.stem

    y = df[args.label_column].astype(int)
    X = df.drop(columns=[args.label_column])
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=args.test_size,
        stratify=y,
        random_state=args.seed,
    )

    model = LogisticRegression(max_iter=args.max_iter, solver="lbfgs")
    model.fit(X_train, y_train)

    train_proba = model.predict_proba(X_train)[:, 1]
    train_pred = (train_proba >= args.threshold).astype(int)

    test_proba = model.predict_proba(X_test)[:, 1]
    test_pred = (test_proba >= args.threshold).astype(int)

    metrics = pd.DataFrame([
        {
            "split": "train",
            "accuracy": accuracy_score(y_train, train_pred),
            "precision": precision_score(y_train, train_pred),
            "recall": recall_score(y_train, train_pred),
            "pr_auc": average_precision_score(y_train, train_proba),
            "roc_auc": roc_auc_score(y_train, train_proba),
        },
        {
            "split": "test",
            "accuracy": accuracy_score(y_test, test_pred),
            "precision": precision_score(y_test, test_pred),
            "recall": recall_score(y_test, test_pred),
            "pr_auc": average_precision_score(y_test, test_proba),
            "roc_auc": roc_auc_score(y_test, test_proba),
        },
    ])

    metrics.to_csv(outdir / f"{version}_metrics.tsv", sep="\t", index=False)

    cm = confusion_matrix(y_test, test_pred)
    pd.DataFrame(
        cm,
        index=["true_0", "true_1"],
        columns=["pred_0", "pred_1"],
    ).to_csv(outdir / f"{version}_confusion_matrix.tsv", sep="\t")

    pd.DataFrame({
        "y_true": y_test.to_numpy(),
        "y_proba": test_proba,
        "y_pred": test_pred,
    }).to_csv(outdir / f"{version}_predictions.csv", index=False)

    weights = pd.DataFrame({
        "feature": X.columns,
        "beta": model.coef_[0],
        "odds_ratio": np.exp(model.coef_[0]),
    }).sort_values("beta", ascending=False)

    weights.to_csv(outdir / f"{version}_feature_weights.csv", index=False)

    if args.permutation:
        perm = permutation_importance(
            model,
            X_test,
            y_test,
            scoring="average_precision",
            n_repeats=args.permutation_repeats,
            random_state=args.seed,
            n_jobs=-1,
        )

        pd.DataFrame({
            "feature": X.columns,
            "importance_mean": perm.importances_mean,
            "importance_std": perm.importances_std,
        }).sort_values("importance_mean", ascending=False).to_csv(
            outdir / f"{version}_permutation_importance.csv",
            index=False,
        )

    print(f"Saved evaluation outputs to: {outdir}")
    print(metrics.to_string(index=False))


def compare_predictions(args):
    import matplotlib.pyplot as plt

    base = pd.read_csv(args.baseline)
    new = pd.read_csv(args.new)

    y_base = base[args.true_column]
    p_base = base[args.proba_column]

    y_new = new[args.true_column]
    p_new = new[args.proba_column]

    precision_base, recall_base, _ = precision_recall_curve(y_base, p_base)
    precision_new, recall_new, _ = precision_recall_curve(y_new, p_new)

    fpr_base, tpr_base, _ = roc_curve(y_base, p_base)
    fpr_new, tpr_new, _ = roc_curve(y_new, p_new)

    pr_auc_base = auc(recall_base, precision_base)
    pr_auc_new = auc(recall_new, precision_new)

    roc_auc_base = auc(fpr_base, tpr_base)
    roc_auc_new = auc(fpr_new, tpr_new)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].plot(recall_base, precision_base, label=f"{args.baseline_label} AUC={pr_auc_base:.3f}")
    axes[0].plot(recall_new, precision_new, label=f"{args.new_label} AUC={pr_auc_new:.3f}")
    axes[0].set_xlabel("Recall")
    axes[0].set_ylabel("Precision")
    axes[0].set_title("Precision recall")
    axes[0].legend()

    axes[1].plot(fpr_base, tpr_base, label=f"{args.baseline_label} AUC={roc_auc_base:.3f}")
    axes[1].plot(fpr_new, tpr_new, label=f"{args.new_label} AUC={roc_auc_new:.3f}")
    axes[1].set_xlabel("False positive rate")
    axes[1].set_ylabel("True positive rate")
    axes[1].set_title("ROC")
    axes[1].legend()

    plt.tight_layout()
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved comparison figure: {args.output}")


def intersect_mpra(args):
    try:
        from pybedtools import BedTool
    except ImportError as exc:
        raise ImportError(
            "pybedtools is required for intersect-mpra. Install it with conda or pip."
        ) from exc

    args.output.parent.mkdir(parents=True, exist_ok=True)

    motifs = BedTool(str(args.motifs))
    regions = BedTool(str(args.regions))

    result = motifs.intersect(regions, wa=True).sort().unique()
    result.saveas(str(args.output))

    print(f"Saved motif-region intersections: {args.output}")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare and evaluate a generic decision table."
    )

    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("aggregate", help="Aggregate cell-specific feature columns.")
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--cells", nargs="+", default=["hepg2", "k562"])
    p.add_argument("--max-tf-binding", type=int, default=3)
    p.set_defaults(func=aggregate_features)

    p = sub.add_parser("align", help="Align a table to a schema/template.")
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--template", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--label", type=int, choices=[0, 1])
    p.add_argument("--label-column", default="activity_score")
    p.set_defaults(func=align_to_schema)

    p = sub.add_parser("balance", help="Balance positive and negative labelled tables.")
    p.add_argument("--positives", required=True, type=Path)
    p.add_argument("--negatives", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--label-column", default="activity_score")
    p.add_argument("--seed", type=int, default=42)
    p.set_defaults(func=balance_table)

    p = sub.add_parser("clean", help="Remove obvious label conflicts.")
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--label-column", default="activity_score")
    p.set_defaults(func=clean_table)

    p = sub.add_parser("evaluate", help="Evaluate a decision table with logistic regression.")
    p.add_argument("--table", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--name", default=None)
    p.add_argument("--label-column", default="activity_score")
    p.add_argument("--test-size", type=float, default=0.2)
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--seed", type=int, default=50)
    p.add_argument("--max-iter", type=int, default=10000)
    p.add_argument("--permutation", action="store_true")
    p.add_argument("--permutation-repeats", type=int, default=30)
    p.set_defaults(func=evaluate_table)

    p = sub.add_parser("compare", help="Compare two prediction files.")
    p.add_argument("--baseline", required=True, type=Path)
    p.add_argument("--new", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--baseline-label", default="Baseline")
    p.add_argument("--new-label", default="New")
    p.add_argument("--true-column", default="y_true")
    p.add_argument("--proba-column", default="y_proba")
    p.add_argument("--dpi", type=int, default=300)
    p.set_defaults(func=compare_predictions)

    p = sub.add_parser("intersect-mpra", help="Intersect motif BED with MPRA region BED.")
    p.add_argument("--motifs", required=True, type=Path)
    p.add_argument("--regions", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.set_defaults(func=intersect_mpra)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
