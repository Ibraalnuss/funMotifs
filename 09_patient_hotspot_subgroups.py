#!/usr/bin/env python3

"""
09_patient_hotspot_subgroups.py

Cluster CRC patients into subgroups based on their recurrent hotspot disruption
profiles, using Jaccard distance and average-linkage hierarchical clustering.

Steps:
  1. Optionally remove hotspots overlapping a genomic blacklist BED.
  2. Filter hotspots to those seen in at least --min-patients patients.
  3. Use bedtools to intersect per-patient motif events with the selected hotspots.
  4. Build a patient x hotspot binary presence/absence matrix.
  5. Compute pairwise Jaccard distances and run average-linkage clustering.
  6. Scan thresholds in a user-defined range; choose the threshold that gives
     ~3 main groups with reasonable sizes.
  7. Assign patients to clusters (C1, C2, ... ordered by size).
  8. Summarize hotspot frequencies, top genes, and top TFs per cluster.
  9. Export gene lists per cluster for downstream pathway enrichment.
 10. Generate figures: dendrogram, heatmap, subgroup profile plots, size barplot.

Requires:
    bedtools in PATH
    scipy for hierarchical clustering

Example:
    python scripts/09_patient_hotspot_subgroups.py \\
      --hotspot-table results/tables/recurrent_hotspots_and_genes/recurrent_hotspot_regions.tsv \\
      --event-table results/tables/variant_overlap_disruption/cleaned/unique_motif_region_per_patient.tsv \\
      --outdir results/tables/patient_hotspot_subgroups \\
      --figdir results/figures/patient_hotspot_subgroups \\
      --min-patients 80
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def require_bedtools() -> None:
    if shutil.which("bedtools") is None:
        raise RuntimeError(
            "bedtools was not found in PATH. Activate the correct conda environment first."
        )


def clean_split(value) -> list[str]:
    if pd.isna(value):
        return []
    value = str(value)
    parts = re.split(r"[;,]", value)
    return [p.strip() for p in parts if p.strip() and p.strip().lower() != "nan"]


def clean_text(value) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def shorten(value, max_len=34) -> str:
    value = clean_text(value)
    if len(value) <= max_len:
        return value
    return value[: max_len - 3] + "..."


def first_gene(value) -> str:
    parts = clean_split(value)
    return parts[0] if parts else ""


def first_tfs(value, max_items=2) -> str:
    parts = clean_split(value)
    if not parts:
        return ""
    if len(parts) <= max_items:
        return ", ".join(parts)
    return ", ".join(parts[:max_items]) + "..."


def parse_thresholds(value: str) -> list[float]:
    value = str(value)

    if ":" in value:
        start, stop, step = [float(x) for x in value.split(":")]
        thresholds = np.arange(start, stop + 1e-9, step)
        return [round(float(x), 4) for x in thresholds]

    return [float(x.strip()) for x in value.split(",") if x.strip()]


def parse_focus_clusters(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    mkdir(path.parent)
    df.to_csv(path, sep="\t", index=False)


def save_text(path: Path, text: str) -> None:
    mkdir(path.parent)
    path.write_text(text)


def read_hotspot_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Hotspot table not found: {path}")

    df = pd.read_csv(path, sep="\t", low_memory=False)
    df["Hotspot_ID"] = df["Hotspot_ID"].astype(str)

    required = [
        "Hotspot_ID",
        "Chr",
        "Hotspot_Start",
        "Hotspot_End",
        "Max_Unique_Patients",
        "Top_TFs",
        "Top_Genes",
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in hotspot table: {missing}")

    return df


def read_event_table(path: Path) -> pd.DataFrame:
    """
    Reads patient-level motif event table and fixes a known merged-header issue.

    Expected columns:
    Chr, Motif_Start, Motif_End, Tumor_Sample_Barcode
    """
    if not path.exists():
        raise FileNotFoundError(f"Event table not found: {path}")

    with path.open("r") as handle:
        header_line = handle.readline().rstrip("\n")
        first_data_line = handle.readline().rstrip("\n")

    header = header_line.split("\t")
    first_data = first_data_line.split("\t")

    if len(header) != len(first_data):
        fixed_header = []

        for col in header:
            if col in {
                "Tumor_Seq_AlleleTumor_Sample_Barcode",
                "Tumor_Seq_Allele2Tumor_Sample_Barcode",
            }:
                fixed_header.extend(["Tumor_Seq_Allele2", "Tumor_Sample_Barcode"])
            else:
                fixed_header.append(col)

        if len(fixed_header) != len(first_data):
            raise ValueError(
                "Could not fix event table header mismatch. "
                f"Header columns: {len(header)}, fixed header columns: {len(fixed_header)}, "
                f"data columns: {len(first_data)}"
            )

        df = pd.read_csv(
            path,
            sep="\t",
            names=fixed_header,
            skiprows=1,
            low_memory=False,
        )
    else:
        df = pd.read_csv(path, sep="\t", low_memory=False)

    required = ["Chr", "Motif_Start", "Motif_End", "Tumor_Sample_Barcode"]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns in event table: {missing}")

    df["Tumor_Sample_Barcode"] = df["Tumor_Sample_Barcode"].astype(str)

    return df


def write_hotspot_bed(df: pd.DataFrame, out_bed: Path) -> None:
    bed = df[["Chr", "Hotspot_Start", "Hotspot_End", "Hotspot_ID"]].copy()
    bed["Hotspot_Start"] = bed["Hotspot_Start"].astype(int) - 1
    bed["Hotspot_End"] = bed["Hotspot_End"].astype(int)

    if (bed["Hotspot_Start"] < 0).any():
        raise ValueError("Negative hotspot BED start detected.")

    bed.to_csv(out_bed, sep="\t", header=False, index=False)


def write_event_bed(df: pd.DataFrame, out_bed: Path) -> None:
    event = df[["Chr", "Motif_Start", "Motif_End", "Tumor_Sample_Barcode"]].copy()
    event = event.dropna()

    event["Motif_Start"] = event["Motif_Start"].astype(int) - 1
    event["Motif_End"] = event["Motif_End"].astype(int)
    event["Tumor_Sample_Barcode"] = event["Tumor_Sample_Barcode"].astype(str)

    if (event["Motif_Start"] < 0).any():
        raise ValueError("Negative motif BED start detected.")

    event.to_csv(out_bed, sep="\t", header=False, index=False)


def run_intersect(event_bed: Path, hotspot_bed: Path, out_file: Path) -> None:
    with out_file.open("w") as out:
        subprocess.run(
            [
                "bedtools",
                "intersect",
                "-a",
                str(event_bed),
                "-b",
                str(hotspot_bed),
                "-wa",
                "-wb",
            ],
            check=True,
            stdout=out,
        )


def filter_blacklist(hotspots: pd.DataFrame, blacklist_bed: Path | None, outdir: Path) -> pd.DataFrame:
    if blacklist_bed is None:
        return hotspots

    if not blacklist_bed.exists():
        raise FileNotFoundError(f"Blacklist BED not found: {blacklist_bed}")

    require_bedtools()

    tmp_dir = outdir / "blacklist_filtering"
    mkdir(tmp_dir)

    input_bed = tmp_dir / "hotspots.for_blacklist_overlap.bed"
    kept_bed = tmp_dir / "hotspots.non_blacklisted.bed"
    removed_bed = tmp_dir / "hotspots.blacklisted_removed.bed"

    write_hotspot_bed(hotspots, input_bed)

    with kept_bed.open("w") as out:
        subprocess.run(
            ["bedtools", "intersect", "-a", str(input_bed), "-b", str(blacklist_bed), "-v"],
            check=True,
            stdout=out,
        )

    with removed_bed.open("w") as out:
        subprocess.run(
            ["bedtools", "intersect", "-a", str(input_bed), "-b", str(blacklist_bed), "-wa"],
            check=True,
            stdout=out,
        )

    kept_ids = read_bed_hotspot_ids(kept_bed)
    removed_ids = read_bed_hotspot_ids(removed_bed)

    kept = hotspots[hotspots["Hotspot_ID"].isin(kept_ids)].copy()
    removed = hotspots[hotspots["Hotspot_ID"].isin(removed_ids)].copy()

    save_tsv(kept, tmp_dir / "hotspots.non_blacklisted.tsv")
    save_tsv(removed, tmp_dir / "hotspots.blacklisted_removed.tsv")

    summary = [
        "Blacklist filtering summary",
        "",
        f"Input hotspots: {len(hotspots):,}",
        f"Kept hotspots: {len(kept):,}",
        f"Removed hotspots: {len(removed):,}",
        f"Blacklist BED: {blacklist_bed}",
    ]
    save_text(tmp_dir / "blacklist_filtering_summary.txt", "\n".join(summary))

    return kept


def read_bed_hotspot_ids(path: Path) -> set[str]:
    if not path.exists() or path.stat().st_size == 0:
        return set()

    bed = pd.read_csv(path, sep="\t", header=None)

    if bed.shape[1] < 4:
        raise ValueError(f"Expected at least 4 columns in BED file: {path}")

    return set(bed.iloc[:, 3].astype(str))


def read_intersections(intersect_file: Path) -> pd.DataFrame:
    if not intersect_file.exists() or intersect_file.stat().st_size == 0:
        raise ValueError("No event-hotspot overlaps found.")

    cols = [
        "Event_Chr",
        "Event_Start_0based",
        "Event_End",
        "Tumor_Sample_Barcode",
        "Hotspot_Chr",
        "Hotspot_Start_0based",
        "Hotspot_End",
        "Hotspot_ID",
    ]

    df = pd.read_csv(intersect_file, sep="\t", header=None, names=cols)
    df["Tumor_Sample_Barcode"] = df["Tumor_Sample_Barcode"].astype(str)
    df["Hotspot_ID"] = df["Hotspot_ID"].astype(str)

    patient_hotspot = (
        df[["Tumor_Sample_Barcode", "Hotspot_ID"]]
        .drop_duplicates()
        .sort_values(["Tumor_Sample_Barcode", "Hotspot_ID"])
        .reset_index(drop=True)
    )

    return patient_hotspot


def make_binary_matrix(patient_hotspot: pd.DataFrame, selected_hotspot_ids: list[str]) -> pd.DataFrame:
    matrix = pd.crosstab(
        patient_hotspot["Tumor_Sample_Barcode"],
        patient_hotspot["Hotspot_ID"],
    )

    matrix = (matrix > 0).astype(int)
    matrix.index.name = "Tumor_Sample_Barcode"

    selected_hotspot_ids = [h for h in selected_hotspot_ids if h in matrix.columns]
    matrix = matrix[selected_hotspot_ids]

    return matrix


def cluster_at_threshold(linkage_matrix, threshold):
    return fcluster(linkage_matrix, t=threshold, criterion="distance")


def summarize_cluster_sizes(cluster_labels):
    return pd.Series(cluster_labels).value_counts().sort_values(ascending=False)


def score_threshold(
    cluster_labels,
    target_main_groups,
    min_main_group_size,
    min_main_group_fraction,
) -> dict:
    sizes = summarize_cluster_sizes(cluster_labels)
    total = sizes.sum()

    main = sizes[sizes >= min_main_group_size]
    n_main = len(main)
    main_fraction = main.sum() / total if total > 0 else 0

    n_clusters = len(sizes)
    n_singletons = int((sizes == 1).sum())
    n_tiny = int((sizes < min_main_group_size).sum())

    largest_fraction = sizes.iloc[0] / total
    top3_fraction = sizes.iloc[:3].sum() / total if len(sizes) >= 3 else sizes.sum() / total

    score = 0.0
    score += 100.0 - 25.0 * abs(n_main - target_main_groups)
    score += 50.0 * main_fraction
    score += 20.0 * top3_fraction
    score -= 2.0 * n_tiny
    score -= 1.0 * n_singletons

    if largest_fraction > 0.90:
        score -= 30.0

    if main_fraction < min_main_group_fraction:
        score -= 40.0

    return {
        "Score": score,
        "Num_Clusters": n_clusters,
        "Num_Main_Groups": n_main,
        "Main_Group_Fraction": main_fraction,
        "Largest_Group_Size": int(sizes.iloc[0]),
        "Largest_Group_Fraction": largest_fraction,
        "Top3_Fraction": top3_fraction,
        "Num_Tiny_Clusters": n_tiny,
        "Num_Singletons": n_singletons,
        "Cluster_Sizes": ";".join(map(str, sizes.tolist())),
    }


def scan_thresholds(
    linkage_matrix,
    thresholds,
    target_main_groups,
    min_main_group_size,
    min_main_group_fraction,
) -> pd.DataFrame:
    records = []

    for threshold in thresholds:
        labels = cluster_at_threshold(linkage_matrix, threshold)
        stats = score_threshold(
            labels,
            target_main_groups,
            min_main_group_size,
            min_main_group_fraction,
        )
        stats["Jaccard_Threshold"] = threshold
        records.append(stats)

    df = pd.DataFrame(records)
    df = df.sort_values(["Score", "Jaccard_Threshold"], ascending=[False, True])
    return df


def relabel_clusters_by_size(raw_labels):
    s = pd.Series(raw_labels)
    size_order = s.value_counts().sort_values(ascending=False).index.tolist()

    mapping = {raw: f"C{i + 1}" for i, raw in enumerate(size_order)}
    new_labels = s.map(mapping).values

    return new_labels, mapping


def build_hotspot_metadata(hotspots: pd.DataFrame) -> dict:
    meta = {}

    for _, row in hotspots.iterrows():
        hs = str(row["Hotspot_ID"])
        meta[hs] = {
            "Hotspot_ID": hs,
            "Top_Genes": row.get("Top_Genes", ""),
            "Top_TFs": row.get("Top_TFs", ""),
            "Representative_Gene": first_gene(row.get("Top_Genes", "")),
            "Representative_TF": first_tfs(row.get("Top_TFs", ""), max_items=2),
            "All_Genes": clean_split(row.get("Top_Genes", "")),
            "All_TFs": clean_split(row.get("Top_TFs", "")),
            "Max_Unique_Patients": row.get("Max_Unique_Patients", np.nan),
            "Chr": row.get("Chr", ""),
            "Hotspot_Start": row.get("Hotspot_Start", ""),
            "Hotspot_End": row.get("Hotspot_End", ""),
        }

    return meta


def compute_cluster_hotspot_profiles(
    matrix: pd.DataFrame,
    cluster_df: pd.DataFrame,
    hotspot_meta: dict,
) -> pd.DataFrame:
    records = []

    for cluster, sub in cluster_df.groupby("Cluster_Label"):
        patients = sub["Tumor_Sample_Barcode"].tolist()
        cluster_matrix = matrix.loc[matrix.index.intersection(patients)]
        size = cluster_matrix.shape[0]

        for hs in matrix.columns:
            count = int(cluster_matrix[hs].sum())
            freq = count / size if size > 0 else 0.0
            meta = hotspot_meta.get(hs, {})

            records.append(
                {
                    "Cluster_Label": cluster,
                    "Cluster_Size": size,
                    "Hotspot_ID": hs,
                    "Patients_With_Hotspot": count,
                    "Frequency": freq,
                    "Percent": 100 * freq,
                    "Representative_Gene": meta.get("Representative_Gene", ""),
                    "Representative_TF": meta.get("Representative_TF", ""),
                    "Top_Genes": meta.get("Top_Genes", ""),
                    "Top_TFs": meta.get("Top_TFs", ""),
                    "Max_Unique_Patients": meta.get("Max_Unique_Patients", np.nan),
                }
            )

    return pd.DataFrame(records)


def summarize_genes_by_cluster(
    profile_df: pd.DataFrame,
    hotspot_meta: dict,
    main_clusters: list[str],
    signature_threshold: float,
) -> pd.DataFrame:
    records = []

    for cl in main_clusters:
        sub_all = profile_df[
            (profile_df["Cluster_Label"] == cl)
            & (profile_df["Frequency"] > 0)
        ].copy()

        sub_sig = profile_df[
            (profile_df["Cluster_Label"] == cl)
            & (profile_df["Frequency"] >= signature_threshold)
        ].copy()

        all_genes = []
        sig_genes = []
        all_tfs = []
        sig_tfs = []

        for hs in sub_all["Hotspot_ID"]:
            all_genes.extend(hotspot_meta.get(hs, {}).get("All_Genes", []))
            all_tfs.extend(hotspot_meta.get(hs, {}).get("All_TFs", []))

        for hs in sub_sig["Hotspot_ID"]:
            sig_genes.extend(hotspot_meta.get(hs, {}).get("All_Genes", []))
            sig_tfs.extend(hotspot_meta.get(hs, {}).get("All_TFs", []))

        records.append(
            {
                "Cluster_Label": cl,
                "Cluster_Size": int(sub_all["Cluster_Size"].iloc[0]) if not sub_all.empty else 0,
                "Num_All_Present_Genes": len(sorted(set(all_genes))),
                "All_Present_Genes": "; ".join(sorted(set(all_genes))),
                "Num_Signature_Genes": len(sorted(set(sig_genes))),
                "Signature_Genes": "; ".join(sorted(set(sig_genes))),
                "Num_All_Present_TFs": len(sorted(set(all_tfs))),
                "All_Present_TFs": "; ".join(sorted(set(all_tfs))),
                "Num_Signature_TFs": len(sorted(set(sig_tfs))),
                "Signature_TFs": "; ".join(sorted(set(sig_tfs))),
                "Signature_Threshold": signature_threshold,
            }
        )

    return pd.DataFrame(records)


def export_gene_lists(gene_summary: pd.DataFrame, outdir: Path) -> None:
    base = outdir / "gene_lists_for_pathway_enrichment"
    signature_dir = base / "signature_genes"
    all_dir = base / "all_present_genes"

    mkdir(signature_dir)
    mkdir(all_dir)

    for _, row in gene_summary.iterrows():
        cl = row["Cluster_Label"]

        sig_genes = clean_split(row.get("Signature_Genes", ""))
        all_genes = clean_split(row.get("All_Present_Genes", ""))

        (signature_dir / f"{cl}_signature_genes.txt").write_text(
            "\n".join(sig_genes) + ("\n" if sig_genes else "")
        )

        (all_dir / f"{cl}_all_present_genes.txt").write_text(
            "\n".join(all_genes) + ("\n" if all_genes else "")
        )


def choose_main_clusters(cluster_df: pd.DataFrame, top_main_clusters: int, min_main_group_size: int) -> list[str]:
    sizes = (
        cluster_df.groupby("Cluster_Label")["Tumor_Sample_Barcode"]
        .nunique()
        .sort_values(ascending=False)
    )

    sizes = sizes[sizes >= min_main_group_size]
    return sizes.head(top_main_clusters).index.tolist()


def plot_threshold_scan(scan_df: pd.DataFrame, chosen_threshold: float, out_png: Path, dpi: int) -> None:
    fig, ax1 = plt.subplots(figsize=(10, 5))

    ax1.plot(scan_df["Jaccard_Threshold"], scan_df["Score"], marker="o")
    ax1.axvline(chosen_threshold, linestyle="--", color="black", linewidth=1)
    ax1.set_xlabel("Jaccard distance threshold")
    ax1.set_ylabel("Threshold score")
    ax1.set_title("Jaccard threshold scan")

    ax2 = ax1.twinx()
    ax2.plot(
        scan_df["Jaccard_Threshold"],
        scan_df["Num_Main_Groups"],
        marker="s",
        alpha=0.7,
    )
    ax2.set_ylabel("Number of main groups")

    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def get_cluster_sizes(cluster_df: pd.DataFrame) -> pd.Series:
    return (
        cluster_df[["Cluster_Label", "Tumor_Sample_Barcode"]]
        .drop_duplicates()
        .groupby("Cluster_Label")["Tumor_Sample_Barcode"]
        .nunique()
        .sort_values(ascending=False)
    )


def cluster_color_map(cluster_order: list[str]) -> dict:
    palette = [
        "#4C78A8",
        "#F58518",
        "#54A24B",
        "#E45756",
        "#72B7B2",
        "#B279A2",
        "#FF9DA6",
        "#9D755D",
        "#BAB0AC",
        "#A0CBE8",
        "#FFBE7D",
        "#8CD17D",
    ]

    return {cl: palette[i % len(palette)] for i, cl in enumerate(cluster_order)}


def plot_cluster_size_bar(cluster_df: pd.DataFrame, out_png: Path, dpi: int) -> None:
    sizes = get_cluster_sizes(cluster_df)
    order = sizes.index.tolist()
    colors = cluster_color_map(order)
    total = int(sizes.sum())

    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    x = np.arange(len(order))
    bars = ax.bar(x, sizes.values, color=[colors[c] for c in order], edgecolor="black", linewidth=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(order, fontsize=11, fontweight="bold")
    ax.set_ylabel("Patients")
    ax.set_title("Sizes of hotspot-defined patient groups")

    for bar, n in zip(bars, sizes.values):
        pct = 100 * n / total
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(sizes.values) * 0.015,
            f"{n}\n({pct:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_cluster_heatmap(
    matrix: pd.DataFrame,
    cluster_df: pd.DataFrame,
    hotspot_meta: dict,
    main_clusters: list[str],
    out_png: Path,
    dpi: int,
) -> None:
    cluster_order = main_clusters
    records = []

    for cl in cluster_order:
        patients = cluster_df.loc[cluster_df["Cluster_Label"] == cl, "Tumor_Sample_Barcode"]
        sub = matrix.loc[matrix.index.intersection(patients)]
        record = {"Cluster_Label": cl, "Cluster_Size": sub.shape[0]}

        for hs in matrix.columns:
            record[hs] = float(sub[hs].mean()) if sub.shape[0] > 0 else 0.0

        records.append(record)

    freq = pd.DataFrame(records)

    hotspot_cols = [c for c in freq.columns if c not in ["Cluster_Label", "Cluster_Size"]]
    hotspot_cols = (
        freq[hotspot_cols]
        .mean(axis=0)
        .sort_values(ascending=False)
        .index.tolist()
    )

    data = freq[hotspot_cols].T.values
    cluster_labels = freq["Cluster_Label"].tolist()

    fig = plt.figure(figsize=(13.5, 7.2))
    gs = GridSpec(
        nrows=1,
        ncols=3,
        width_ratios=[4.2, 7.0, 0.35],
        wspace=0.05,
        figure=fig,
    )

    ax_labels = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])

    im = ax.imshow(data, aspect="auto", cmap="Blues", vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(cluster_labels)))
    ax.set_xticklabels(cluster_labels, fontsize=11, fontweight="bold")
    ax.set_yticks(np.arange(len(hotspot_cols)))
    ax.set_yticklabels([""] * len(hotspot_cols))
    ax.set_xlabel("Patient groups")
    ax.set_title("Hotspot frequencies across patient groups")

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            text_color = "white" if val >= 0.55 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8.5, color=text_color)

    draw_hotspot_metadata_axis(ax_labels, hotspot_cols, hotspot_meta)

    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Within-cluster frequency")

    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    save_tsv(freq, out_png.with_suffix(".frequency_table.tsv"))


def draw_hotspot_metadata_axis(ax, hotspot_order: list[str], hotspot_meta: dict) -> None:
    n = len(hotspot_order)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, n - 0.5)
    ax.axis("off")

    x_hotspot = 0.02
    x_gene = 0.34
    x_tf = 0.57

    ax.text(x_hotspot, n - 0.05, "Hotspot", fontsize=10, fontweight="bold", ha="left", va="bottom")
    ax.text(x_gene, n - 0.05, "Gene", fontsize=10, fontweight="bold", color="#1F77B4", ha="left", va="bottom")
    ax.text(x_tf, n - 0.05, "TF", fontsize=10, fontweight="bold", color="#8B1E3F", ha="left", va="bottom")

    for row_index, hs in enumerate(reversed(hotspot_order)):
        y = row_index
        meta = hotspot_meta.get(hs, {})

        ax.text(x_hotspot, y, hs, fontsize=9, fontweight="bold", ha="left", va="center")
        ax.text(
            x_gene,
            y,
            shorten(meta.get("Representative_Gene", ""), 18),
            fontsize=9,
            color="#1F77B4",
            ha="left",
            va="center",
        )
        ax.text(
            x_tf,
            y,
            shorten(meta.get("Representative_TF", ""), 28),
            fontsize=9,
            color="#8B1E3F",
            ha="left",
            va="center",
        )


def plot_dendrogram_with_group_strip(
    matrix: pd.DataFrame,
    cluster_df: pd.DataFrame,
    linkage_matrix,
    threshold: float,
    main_clusters: list[str],
    out_png: Path,
    dpi: int,
) -> None:
    patient_order = list(matrix.index)
    n_leaves = len(patient_order)

    cluster_lookup = dict(zip(cluster_df["Tumor_Sample_Barcode"], cluster_df["Cluster_Label"]))
    color_map = cluster_color_map(main_clusters + ["Other"])
    color_map["Other"] = "#BDBDBD"

    leaf_group = {}
    for patient in patient_order:
        cl = cluster_lookup.get(patient, "Other")
        if cl not in main_clusters:
            cl = "Other"
        leaf_group[patient] = cl

    node_to_leaves = {i: [patient_order[i]] for i in range(n_leaves)}
    node_group = {}

    for i, row in enumerate(linkage_matrix):
        left = int(row[0])
        right = int(row[1])
        node_id = n_leaves + i

        leaves = node_to_leaves[left] + node_to_leaves[right]
        node_to_leaves[node_id] = leaves

        groups = {leaf_group[p] for p in leaves}
        node_group[node_id] = list(groups)[0] if len(groups) == 1 else "Mixed"

    def link_color_func(k):
        if k < n_leaves:
            return "#7A7A7A"
        group = node_group.get(k, "Mixed")
        if group in color_map:
            return color_map[group]
        return "#8A8A8A"

    fig = plt.figure(figsize=(16.5, 7.8))
    gs = GridSpec(
        nrows=3,
        ncols=1,
        height_ratios=[5.2, 0.40, 1.15],
        hspace=0.08,
        figure=fig,
    )

    ax_dendro = fig.add_subplot(gs[0, 0])
    ax_strip = fig.add_subplot(gs[1, 0])
    ax_legend = fig.add_subplot(gs[2, 0])

    dendro = dendrogram(
        linkage_matrix,
        no_labels=True,
        color_threshold=threshold,
        above_threshold_color="#2F2F2F",
        link_color_func=link_color_func,
        ax=ax_dendro,
    )

    ax_dendro.axhline(threshold, color="#D62728", linestyle="--", linewidth=1.5)
    ax_dendro.set_title("Patient clustering based on Jaccard distance")
    ax_dendro.set_ylabel("Jaccard distance")
    ax_dendro.set_xticks([])
    ax_dendro.spines["top"].set_visible(False)
    ax_dendro.spines["right"].set_visible(False)

    leaf_indices = dendro["leaves"]
    ordered_patients = [patient_order[i] for i in leaf_indices]
    strip_labels = []

    for patient in ordered_patients:
        cl = cluster_lookup.get(patient, "Other")
        if cl not in main_clusters:
            cl = "Other"
        strip_labels.append(cl)

    label_order = main_clusters + ["Other"]
    label_to_int = {lab: i for i, lab in enumerate(label_order)}
    strip_values = np.array([[label_to_int[x] for x in strip_labels]])

    cmap = plt.matplotlib.colors.ListedColormap([color_map[x] for x in label_order])
    ax_strip.imshow(strip_values, aspect="auto", cmap=cmap, interpolation="nearest")
    ax_strip.set_yticks([])
    ax_strip.set_xticks([])
    ax_strip.set_ylabel("Group", rotation=0, labelpad=30, va="center")

    ax_legend.axis("off")
    handles = [Patch(facecolor=color_map[x], edgecolor="black", label=x) for x in label_order]
    handles.append(Patch(facecolor="#D62728", edgecolor="#D62728", label="Cut threshold"))
    ax_legend.legend(handles=handles, loc="center", ncol=min(5, len(handles)), frameon=False)

    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_clear_subgroup_profiles(
    profile_df: pd.DataFrame,
    main_clusters: list[str],
    out_png: Path,
    dpi: int,
    top_hotspots_per_group: int,
    signature_threshold: float,
) -> None:
    if profile_df.empty:
        return

    focus_tables = {}

    for cl in main_clusters:
        sub = profile_df[profile_df["Cluster_Label"] == cl].copy()
        sub = sub.sort_values(["Frequency", "Max_Unique_Patients"], ascending=[False, False])
        sub = sub.head(top_hotspots_per_group).copy()
        focus_tables[cl] = sub

    n_clusters = len(main_clusters)
    fig_height = max(5.5, 3.0 * n_clusters)
    fig = plt.figure(figsize=(15, fig_height))
    gs = GridSpec(nrows=n_clusters, ncols=2, width_ratios=[4.2, 6.0], hspace=0.65, figure=fig)

    for i, cl in enumerate(main_clusters):
        table = focus_tables[cl]
        ax_labels = fig.add_subplot(gs[i, 0])
        ax = fig.add_subplot(gs[i, 1])

        draw_profile_label_columns(ax_labels, table, show_header=(i == 0))

        y = np.arange(len(table))
        ax.barh(y, table["Percent"], color="#4C78A8", edgecolor="black", linewidth=0.4)
        ax.axvline(signature_threshold * 100, linestyle="--", color="#D62728", linewidth=1.2)
        ax.set_yticks(y)
        ax.set_yticklabels([""] * len(y))
        ax.invert_yaxis()
        ax.set_xlim(0, 100)
        ax.set_xlabel("Patients in subgroup with hotspot (%)")
        ax.set_title(f"{cl} hotspot profile, n={int(table['Cluster_Size'].iloc[0]) if not table.empty else 0}")

        for yi, (_, row) in enumerate(table.iterrows()):
            ax.text(
                row["Percent"] + 1,
                yi,
                f"{row['Percent']:.0f}% ({int(row['Patients_With_Hotspot'])}/{int(row['Cluster_Size'])})",
                va="center",
                fontsize=9,
            )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def draw_profile_label_columns(ax, table: pd.DataFrame, show_header: bool = False) -> None:
    n = len(table)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()
    ax.axis("off")

    x_hs = 0.02
    x_gene = 0.36
    x_tf = 0.58

    if show_header:
        ax.text(x_hs, 1.04, "Hotspot", fontsize=10.5, fontweight="bold", transform=ax.transAxes)
        ax.text(x_gene, 1.04, "Gene", fontsize=10.5, fontweight="bold", color="#1F77B4", transform=ax.transAxes)
        ax.text(x_tf, 1.04, "TF", fontsize=10.5, fontweight="bold", color="#8B1E3F", transform=ax.transAxes)

    for y, (_, row) in enumerate(table.iterrows()):
        ax.text(x_hs, y, row["Hotspot_ID"], fontsize=9.2, fontweight="bold", ha="left", va="center")
        ax.text(
            x_gene,
            y,
            shorten(row["Representative_Gene"], 16),
            fontsize=9.2,
            color="#1F77B4",
            ha="left",
            va="center",
        )
        ax.text(
            x_tf,
            y,
            shorten(row["Representative_TF"], 26),
            fontsize=9.2,
            color="#8B1E3F",
            ha="left",
            va="center",
        )


def run(args: argparse.Namespace) -> None:
    require_bedtools()

    mkdir(args.outdir)
    mkdir(args.figdir)

    print("Loading hotspot table...")
    hotspots = read_hotspot_table(args.hotspot_table)
    print(f"Input hotspots: {len(hotspots):,}")

    hotspots = filter_blacklist(hotspots, args.blacklist_bed, args.outdir)

    hotspots["Max_Unique_Patients"] = pd.to_numeric(
        hotspots["Max_Unique_Patients"],
        errors="coerce",
    )

    selected = hotspots[hotspots["Max_Unique_Patients"] >= args.min_patients].copy()
    selected = selected.sort_values(
        ["Max_Unique_Patients", "Hotspot_ID"],
        ascending=[False, True],
    ).head(args.max_hotspots).copy()

    if selected.empty:
        raise RuntimeError("No hotspots remained after filtering.")

    selected_ids = selected["Hotspot_ID"].astype(str).tolist()

    selected_file = args.outdir / f"selected_hotspots_ge{args.min_patients}.tsv"
    save_tsv(selected, selected_file)

    print(f"Selected hotspots: {len(selected):,}")
    print(f"Saved selected hotspots: {selected_file}")

    print("Loading event table...")
    events = read_event_table(args.event_table)
    print(f"Input events: {len(events):,}")

    hotspot_bed = args.outdir / "selected_hotspots.bed"
    event_bed = args.outdir / "patient_motif_events.bed"
    intersect_file = args.outdir / "patient_hotspot_intersections.tsv"

    write_hotspot_bed(selected, hotspot_bed)
    write_event_bed(events, event_bed)
    run_intersect(event_bed, hotspot_bed, intersect_file)

    patient_hotspot = read_intersections(intersect_file)
    save_tsv(patient_hotspot, args.outdir / "patient_hotspot_pairs.tsv")

    matrix = make_binary_matrix(patient_hotspot, selected_ids)

    if matrix.shape[0] < 2:
        raise RuntimeError("Need at least two patients for clustering.")

    if matrix.shape[1] < 1:
        raise RuntimeError("No hotspot columns remained in the binary matrix.")

    save_tsv(matrix.reset_index(), args.outdir / "patient_by_hotspot_binary_matrix.tsv")

    print(f"Matrix shape: {matrix.shape[0]:,} patients x {matrix.shape[1]:,} hotspots")

    print("Computing Jaccard distances and linkage...")
    distances = pdist(matrix.values.astype(bool), metric="jaccard")
    linkage_matrix = linkage(distances, method="average")

    thresholds = parse_thresholds(args.thresholds)
    scan_df = scan_thresholds(
        linkage_matrix=linkage_matrix,
        thresholds=thresholds,
        target_main_groups=args.target_main_groups,
        min_main_group_size=args.min_main_group_size,
        min_main_group_fraction=args.min_main_group_fraction,
    )

    save_tsv(scan_df, args.outdir / "jaccard_threshold_scan.tsv")

    chosen_threshold = float(scan_df.iloc[0]["Jaccard_Threshold"])
    raw_labels = cluster_at_threshold(linkage_matrix, chosen_threshold)
    labels, raw_to_label = relabel_clusters_by_size(raw_labels)

    cluster_df = pd.DataFrame(
        {
            "Tumor_Sample_Barcode": matrix.index.astype(str),
            "Raw_Jaccard_Cluster": raw_labels,
            "Cluster_Label": labels,
            "Jaccard_Threshold": chosen_threshold,
        }
    )

    cluster_df["Jaccard_Cluster"] = cluster_df["Cluster_Label"].str.replace("C", "", regex=False)

    save_tsv(cluster_df, args.outdir / "patient_jaccard_clusters.tsv")

    cluster_sizes = get_cluster_sizes(cluster_df)
    save_tsv(
        cluster_sizes.reset_index().rename(
            columns={"index": "Cluster_Label", "Tumor_Sample_Barcode": "Num_Patients"}
        ),
        args.outdir / "patient_cluster_sizes.tsv",
    )

    hotspot_meta = build_hotspot_metadata(selected)

    profile_df = compute_cluster_hotspot_profiles(matrix, cluster_df, hotspot_meta)
    save_tsv(profile_df, args.outdir / "cluster_hotspot_profiles.tsv")

    main_clusters = choose_main_clusters(
        cluster_df,
        top_main_clusters=args.top_main_clusters,
        min_main_group_size=args.min_main_group_size,
    )

    if args.focus_clusters:
        requested = parse_focus_clusters(args.focus_clusters)
        main_clusters = [cl for cl in requested if cl in set(cluster_df["Cluster_Label"])]

    if not main_clusters:
        main_clusters = cluster_sizes.head(args.top_main_clusters).index.tolist()

    gene_summary = summarize_genes_by_cluster(
        profile_df,
        hotspot_meta,
        main_clusters,
        args.signature_threshold,
    )

    save_tsv(gene_summary, args.outdir / "cluster_gene_tf_summary.tsv")
    export_gene_lists(gene_summary, args.outdir)

    plot_threshold_scan(
        scan_df,
        chosen_threshold,
        args.figdir / "jaccard_threshold_scan.png",
        args.dpi,
    )

    plot_cluster_size_bar(
        cluster_df,
        args.figdir / "patient_group_sizes.png",
        args.dpi,
    )

    plot_cluster_heatmap(
        matrix,
        cluster_df,
        hotspot_meta,
        main_clusters,
        args.figdir / "cluster_hotspot_frequency_heatmap.png",
        args.dpi,
    )

    plot_dendrogram_with_group_strip(
        matrix,
        cluster_df,
        linkage_matrix,
        chosen_threshold,
        main_clusters,
        args.figdir / "patient_jaccard_dendrogram.png",
        args.dpi,
    )

    plot_clear_subgroup_profiles(
        profile_df,
        main_clusters,
        args.figdir / "clear_subgroup_hotspot_profiles.png",
        args.dpi,
        args.top_hotspots_per_group,
        args.signature_threshold,
    )

    summary_lines = [
        "Patient hotspot subgroup summary",
        "",
        f"Input hotspot table: {args.hotspot_table}",
        f"Input event table: {args.event_table}",
        f"Minimum hotspot recurrence: {args.min_patients}",
        f"Maximum selected hotspots: {args.max_hotspots}",
        f"Selected hotspots: {len(selected):,}",
        f"Patients in matrix: {matrix.shape[0]:,}",
        f"Hotspots in matrix: {matrix.shape[1]:,}",
        f"Chosen Jaccard threshold: {chosen_threshold}",
        f"Main clusters: {', '.join(main_clusters)}",
        "",
        "Cluster sizes:",
        cluster_sizes.to_string(),
    ]

    save_text(args.outdir / "patient_hotspot_subgroups_summary.txt", "\n".join(summary_lines))

    print("\nDone.")
    print(f"Tables saved in: {args.outdir}")
    print(f"Figures saved in: {args.figdir}")
    print(f"Chosen Jaccard threshold: {chosen_threshold}")
    print("Main clusters:", ", ".join(main_clusters))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Discover patient subgroups from recurrent disrupted hotspot profiles."
    )

    parser.add_argument(
        "--hotspot-table",
        required=True,
        type=Path,
        help="Recurrent hotspot table from script 06."
    )
    parser.add_argument(
        "--event-table",
        required=True,
        type=Path,
        help="Patient-level motif event table, usually unique_motif_region_per_patient.tsv."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for subgroup tables."
    )
    parser.add_argument(
        "--figdir",
        required=True,
        type=Path,
        help="Output directory for subgroup figures."
    )
    parser.add_argument(
        "--blacklist-bed",
        type=Path,
        default=None,
        help="Optional hg38 blacklist BED. Hotspots overlapping this file are removed."
    )
    parser.add_argument(
        "--min-patients",
        type=int,
        default=80,
        help="Minimum Max_Unique_Patients required for hotspot selection."
    )
    parser.add_argument(
        "--max-hotspots",
        type=int,
        default=50,
        help="Maximum number of selected hotspots after recurrence filtering."
    )
    parser.add_argument(
        "--thresholds",
        default="0.40:0.90:0.02",
        help="Thresholds to scan. Use comma list or start:stop:step."
    )
    parser.add_argument(
        "--target-main-groups",
        type=int,
        default=3,
        help="Preferred number of main patient groups."
    )
    parser.add_argument(
        "--min-main-group-size",
        type=int,
        default=10,
        help="Minimum size for a group to count as a main group."
    )
    parser.add_argument(
        "--min-main-group-fraction",
        type=float,
        default=0.80,
        help="Preferred fraction of patients covered by main groups."
    )
    parser.add_argument(
        "--top-main-clusters",
        type=int,
        default=3,
        help="Number of largest main clusters to summarize."
    )
    parser.add_argument(
        "--focus-clusters",
        default=None,
        help="Optional comma-separated cluster labels to focus on, for example C7,C1,C2."
    )
    parser.add_argument(
        "--top-hotspots-per-group",
        type=int,
        default=5,
        help="Number of hotspots shown per group in profile plots."
    )
    parser.add_argument(
        "--signature-threshold",
        type=float,
        default=0.50,
        help="Minimum within-cluster hotspot frequency for signature genes."
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure resolution."
    )

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
