#!/usr/bin/env python3

"""
08_make_thesis_figures.py

Generate thesis figures from the final analysis tables.

All input paths are passed as command-line arguments -- nothing is hard-coded.
Use --group to select which figures to generate, or pass 'all' to run everything.

Figure groups:
    scoring     Colon motif scoring summary; normalized TF disruption barplot
    hotspots    Hotspot landscape and recurrent loci Manhattan-style plots
    loci        Top integrated recurrent loci barplot
    gwas        GWAS hotspot dotplot and direct-overlap schematic
    tf_gene     TF prioritization barplot and TF-gene bubble plot
    pathways    Pathway enrichment dotplots for top genes and NCG gene sets
    all         Run all groups above

Example:
    python scripts/08_make_thesis_figures.py \\
      --group all \\
      --recurrent-dir results/tables/recurrent_hotspots_and_genes \\
      --external-dir results/tables/external_interpretation \\
      --variant-dir results/tables/variant_overlap_disruption \\
      --outdir results/figures
"""

from __future__ import annotations

import argparse
import math
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


CHR_LENGTHS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}

CHR_ORDER = list(CHR_LENGTHS.keys())


def clean_text(value) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def normalize_chr(chrom) -> str:
    chrom = clean_text(chrom)
    if not chrom:
        return ""
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return chrom


def motif_to_tf(motif_name: str) -> str:
    motif_name = clean_text(motif_name)
    if not motif_name:
        return ""
    if "_" in motif_name:
        return motif_name.rsplit("_", 1)[0].strip()
    return motif_name


def shorten_text(value, max_len=35) -> str:
    value = clean_text(value)
    if len(value) <= max_len:
        return value
    return value[: max_len - 3] + "..."


def split_items(text) -> list[str]:
    text = clean_text(text)
    if not text:
        return []
    return [x.strip() for x in text.split(",") if x.strip()]


def scale_sizes(values, min_size=120, max_size=1400):
    values = pd.Series(values).astype(float)

    if values.empty:
        return values

    vmin = values.min()
    vmax = values.max()

    if vmin == vmax:
        return pd.Series([0.5 * (min_size + max_size)] * len(values), index=values.index)

    return min_size + (values - vmin) * (max_size - min_size) / (vmax - vmin)


def make_offsets():
    offsets = {}
    tick_positions = []
    tick_labels = []

    running = 0

    for chrom in CHR_ORDER:
        offsets[chrom] = running
        tick_positions.append(running + CHR_LENGTHS[chrom] / 2)
        tick_labels.append(chrom.replace("chr", ""))
        running += CHR_LENGTHS[chrom]

    return offsets, tick_positions, tick_labels


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def save_fig(fig, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def maybe_read(path: Path | None, sep="\t") -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep=sep, low_memory=False)


def check_columns(df: pd.DataFrame, required: list[str], label: str) -> bool:
    missing = [col for col in required if col not in df.columns]
    if missing:
        print(f"[SKIP] {label}: missing columns {missing}")
        return False
    return True


def scoring_summary_figure(total: int, functional: int, outdir: Path, dpi: int) -> None:
    if total <= 0 or functional < 0:
        print("[SKIP] scoring summary: invalid counts")
        return

    nonfunctional = total - functional

    df = pd.DataFrame(
        {
            "Category": [
                "Total motifs scored",
                "Predicted functional",
                "Predicted non-functional",
            ],
            "Count": [total, functional, nonfunctional],
        }
    )
    df["Percent_of_total"] = df["Count"] / total * 100

    save_tsv(df, outdir / "section5_scoring_summary_table.tsv")

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(df["Category"], df["Count"])
    ax.set_ylabel("Number of motif instances")
    ax.set_title("Colon motif scoring summary")

    for i, row in df.iterrows():
        ax.text(
            i,
            row["Count"] + total * 0.01,
            f"{int(row['Count']):,}\n({row['Percent_of_total']:.1f}%)",
            ha="center",
            va="bottom",
        )

    ax.tick_params(axis="x", rotation=10)
    fig.tight_layout()
    save_fig(fig, outdir / "section5_scoring_summary.png", dpi)


def top_normalized_tfs_figure(tf_background_file: Path | None, outdir: Path, dpi: int, top_n: int) -> None:
    df = maybe_read(tf_background_file)

    if df.empty:
        print("[SKIP] normalized TF plot: no input file")
        return

    required = ["TF", "Patients_per_1000_Functional_Motifs"]
    if not check_columns(df, required, "normalized TF plot"):
        return

    plot_df = (
        df.sort_values("Patients_per_1000_Functional_Motifs", ascending=False)
        .head(top_n)
        .sort_values("Patients_per_1000_Functional_Motifs", ascending=True)
        .copy()
    )

    save_tsv(plot_df, outdir / "top_normalized_tfs_table.tsv")

    fig, ax = plt.subplots(figsize=(11, 8))
    ax.barh(plot_df["TF"], plot_df["Patients_per_1000_Functional_Motifs"])
    ax.set_xlabel("Unique patients per 1,000 functional motif instances")
    ax.set_ylabel("TF")
    ax.set_title("Top normalized disrupted TF motifs")

    maxv = plot_df["Patients_per_1000_Functional_Motifs"].max()
    for i, v in enumerate(plot_df["Patients_per_1000_Functional_Motifs"]):
        ax.text(v + maxv * 0.01, i, f"{v:.2f}", va="center")

    fig.tight_layout()
    save_fig(fig, outdir / "top_normalized_tfs_barplot.png", dpi)


def tf_background_comparison(functional_table: Path | None, disruptions_file: Path | None, outdir: Path) -> Path | None:
    fun = maybe_read(functional_table)
    dis = maybe_read(disruptions_file)

    if fun.empty or dis.empty:
        print("[SKIP] TF background comparison: missing input files")
        return None

    if "name" not in fun.columns:
        print("[SKIP] TF background comparison: functional table needs column 'name'")
        return None

    if "Name" not in dis.columns or "Tumor_Sample_Barcode" not in dis.columns:
        print("[SKIP] TF background comparison: disruption table needs Name and Tumor_Sample_Barcode")
        return None

    fun["TF"] = fun["name"].map(motif_to_tf)
    dis["TF"] = dis["Name"].map(motif_to_tf)

    bg = fun.groupby("TF").size().reset_index(name="Functional_Motif_Instances")
    patients = (
        dis.groupby("TF")["Tumor_Sample_Barcode"]
        .nunique()
        .reset_index(name="Unique_Patients_Disrupted")
    )

    out = patients.merge(bg, on="TF", how="left")
    out["Functional_Motif_Instances"] = out["Functional_Motif_Instances"].fillna(0).astype(int)

    out["Patients_per_1000_Functional_Motifs"] = out.apply(
        lambda r: (
            r["Unique_Patients_Disrupted"] / r["Functional_Motif_Instances"] * 1000
            if r["Functional_Motif_Instances"] > 0
            else np.nan
        ),
        axis=1,
    )

    out = out.sort_values("Patients_per_1000_Functional_Motifs", ascending=False)
    path = outdir / "tf_background_comparison.tsv"
    save_tsv(out, path)
    return path


def top_integrated_loci_figure(integrated_file: Path | None, outdir: Path, dpi: int, top_n: int) -> None:
    df = maybe_read(integrated_file)

    if df.empty:
        print("[SKIP] top integrated loci: no input file")
        return

    required = ["Chr", "Motif_Start", "Motif_End", "Unique_Patients_Locus", "Top_TFs", "Top_Genes"]
    if not check_columns(df, required, "top integrated loci"):
        return

    sort_cols = ["Unique_Patients_Locus"]
    asc = [False]

    for col in ["Unique_Variants_Locus", "Num_Disrupting_Rows_Locus"]:
        if col in df.columns:
            sort_cols.append(col)
            asc.append(False)

    plot_df = df.sort_values(sort_cols, ascending=asc).head(top_n).copy()
    plot_df["Locus_Label"] = plot_df.apply(
        lambda r: f"{r['Chr']}:{int(r['Motif_Start']):,}-{int(r['Motif_End']):,}",
        axis=1,
    )
    plot_df["Annotation"] = plot_df.apply(
        lambda r: (
            f"TF: {shorten_text(r.get('Top_TFs', ''), 32) or 'NA'}\n"
            f"Gene: {shorten_text(r.get('Top_Genes', ''), 32) or 'NA'}"
        ),
        axis=1,
    )

    save_tsv(plot_df, outdir / "top_integrated_recurrent_loci_plotted.tsv")

    plot_df = plot_df.iloc[::-1].copy()
    y = range(len(plot_df))
    counts = plot_df["Unique_Patients_Locus"].astype(int)
    max_count = max(counts) if len(counts) else 1

    fig, ax = plt.subplots(figsize=(16, 9))
    ax.barh(y, counts)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["Locus_Label"], fontsize=10)
    ax.set_xlabel("Unique patients")
    ax.set_ylabel("Recurrent disrupted locus")
    ax.set_title("Top integrated recurrent loci")
    ax.set_xlim(0, max_count * 1.75)
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)

    count_x = max_count * 1.05
    annotation_x = max_count * 1.18

    for i, (_, row) in enumerate(plot_df.iterrows()):
        ax.text(count_x, i, str(int(row["Unique_Patients_Locus"])), va="center", fontsize=10)
        ax.text(annotation_x, i, row["Annotation"], va="center", fontsize=8)

    fig.tight_layout()
    save_fig(fig, outdir / "top_integrated_recurrent_loci.png", dpi)


def recurrent_loci_manhattan(
    loci_file: Path | None,
    outdir: Path,
    dpi: int,
    min_patients: int,
    top_labels: int,
) -> None:
    df = maybe_read(loci_file)

    if df.empty:
        print("[SKIP] recurrent loci Manhattan: no input file")
        return

    if "Unique_Patients" not in df.columns and "Unique_Patients_Locus" in df.columns:
        df = df.rename(columns={"Unique_Patients_Locus": "Unique_Patients"})

    required = ["Chr", "Motif_Start", "Motif_End", "Unique_Patients"]
    if not check_columns(df, required, "recurrent loci Manhattan"):
        return

    df["Chr"] = df["Chr"].map(normalize_chr)
    df = df[df["Chr"].isin(CHR_ORDER)].copy()
    df = df[pd.to_numeric(df["Unique_Patients"], errors="coerce") >= min_patients].copy()

    if df.empty:
        print("[SKIP] recurrent loci Manhattan: no rows after filtering")
        return

    df["Locus_Midpoint"] = ((df["Motif_Start"].astype(int) + df["Motif_End"].astype(int)) / 2).astype(int)

    offsets, tick_positions, tick_labels = make_offsets()
    df["Genome_X"] = df.apply(lambda r: offsets[r["Chr"]] + r["Locus_Midpoint"], axis=1)

    save_tsv(df.sort_values(["Chr", "Locus_Midpoint"]), outdir / "recurrent_loci_manhattan_input.tsv")

    chrom_to_color_index = {chrom: i % 2 for i, chrom in enumerate(CHR_ORDER)}
    colors = df["Chr"].map(lambda c: "#4C78A8" if chrom_to_color_index[c] == 0 else "#F58518")

    fig, ax = plt.subplots(figsize=(16, 6))
    ax.scatter(df["Genome_X"], df["Unique_Patients"], c=colors, s=18, alpha=0.75)

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Unique patients")
    ax.set_title("Genome-wide recurrent disrupted motif loci")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)

    labels = df.sort_values("Unique_Patients", ascending=False).head(top_labels)
    for _, row in labels.iterrows():
        label = f"{row['Chr']}:{int(row['Motif_Start'])/1e6:.2f} Mb"
        ax.text(row["Genome_X"], row["Unique_Patients"], label, fontsize=8, rotation=25)

    fig.tight_layout()
    save_fig(fig, outdir / "manhattan_recurrent_loci.png", dpi)


def hotspot_landscape_figure(
    hotspot_file: Path | None,
    outdir: Path,
    dpi: int,
    top_labels: int,
) -> None:
    df = maybe_read(hotspot_file)

    if df.empty:
        print("[SKIP] hotspot landscape: no input file")
        return

    required = ["Hotspot_ID", "Chr", "Hotspot_Start", "Hotspot_End", "Max_Unique_Patients"]
    if not check_columns(df, required, "hotspot landscape"):
        return

    df = df.copy()
    df["Chr"] = df["Chr"].map(normalize_chr)
    df["Chromosome"] = df["Chr"].str.replace("chr", "", regex=False)
    df["Chromosome_Num"] = pd.to_numeric(df["Chromosome"], errors="coerce")

    df["Hotspot_Start"] = pd.to_numeric(df["Hotspot_Start"], errors="coerce")
    df["Hotspot_End"] = pd.to_numeric(df["Hotspot_End"], errors="coerce")
    df["Max_Unique_Patients"] = pd.to_numeric(df["Max_Unique_Patients"], errors="coerce")

    df = df.dropna(subset=["Chromosome_Num", "Hotspot_Start", "Hotspot_End", "Max_Unique_Patients"]).copy()
    df = df[(df["Chromosome_Num"] >= 1) & (df["Chromosome_Num"] <= 22)].copy()

    if df.empty:
        print("[SKIP] hotspot landscape: no autosomal rows")
        return

    df["Midpoint"] = ((df["Hotspot_Start"] + df["Hotspot_End"]) / 2).astype(int)

    if "Hotspot_Width_bp" not in df.columns:
        df["Hotspot_Width_bp"] = df["Hotspot_End"] - df["Hotspot_Start"] + 1
    else:
        df["Hotspot_Width_bp"] = pd.to_numeric(df["Hotspot_Width_bp"], errors="coerce")
        df["Hotspot_Width_bp"] = df["Hotspot_Width_bp"].fillna(df["Hotspot_End"] - df["Hotspot_Start"] + 1)

    chrom_dfs = []
    offset = 0
    xticks = []
    xticklabels = []

    for chrom in range(1, 23):
        sub = df[df["Chromosome_Num"] == chrom].sort_values("Midpoint").copy()
        if sub.empty:
            continue

        span = sub["Midpoint"].max() - sub["Midpoint"].min() + 1
        sub["Plot_X"] = sub["Midpoint"] - sub["Midpoint"].min() + offset
        chrom_dfs.append(sub)
        xticks.append(offset + span / 2)
        xticklabels.append(str(chrom))
        offset += span + 8_000_000

    plot_df = pd.concat(chrom_dfs, ignore_index=True)
    plot_df["Color"] = plot_df["Chromosome_Num"].apply(lambda x: "#8ea8cf" if int(x) % 2 else "#b8b8b8")
    plot_df["Point_Size"] = scale_sizes(np.sqrt(plot_df["Hotspot_Width_bp"].clip(lower=1)), 8, 80)

    save_tsv(plot_df, outdir / "hotspot_landscape_manhattan_input.tsv")

    fig, ax = plt.subplots(figsize=(16, 6))
    ax.scatter(
        plot_df["Plot_X"],
        plot_df["Max_Unique_Patients"],
        s=plot_df["Point_Size"],
        c=plot_df["Color"],
        alpha=0.8,
        edgecolors="none",
    )

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Maximum unique patients")
    ax.set_title("Recurrent disrupted hotspot landscape")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)

    labels = plot_df.sort_values("Max_Unique_Patients", ascending=False).head(top_labels)
    for _, row in labels.iterrows():
        label = f"{row['Chr']}:{int(row['Midpoint'])/1e6:.2f} Mb\n{int(row['Hotspot_Width_bp'])} bp"
        ax.text(row["Plot_X"], row["Max_Unique_Patients"], label, fontsize=8)

    fig.tight_layout()
    save_fig(fig, outdir / "hotspot_landscape_manhattan_with_width.png", dpi)


def gwas_hotspot_dotplot(gwas_besthit_file: Path | None, outdir: Path, dpi: int, top_n: int) -> None:
    df = maybe_read(gwas_besthit_file)

    if df.empty:
        print("[SKIP] GWAS dotplot: no input file")
        return

    required = [
        "Top_Genes",
        "Top_TFs",
        "Max_Unique_Patients",
        "GWAS_SNP",
        "Distance_to_GWAS_bp",
        "Evidence_Tier",
        "GWAS_Trait",
        "Gene_Concordance",
    ]
    if not check_columns(df, required, "GWAS dotplot"):
        return

    tier_order = {
        "Tier 1: exact": 0,
        "Tier 1: <=1kb": 1,
        "Tier 2: <=10kb": 2,
        "Tier 3: <=50kb": 3,
    }

    df = df.copy()
    df["Evidence_Tier_Rank"] = df["Evidence_Tier"].map(tier_order).fillna(99)
    df["Distance_to_GWAS_bp"] = pd.to_numeric(df["Distance_to_GWAS_bp"], errors="coerce")
    df["Max_Unique_Patients"] = pd.to_numeric(df["Max_Unique_Patients"], errors="coerce")

    df = (
        df.sort_values(["Distance_to_GWAS_bp", "Gene_Concordance", "Max_Unique_Patients"], ascending=[True, False, False])
        .head(top_n)
        .copy()
    )

    df["Gene_Label"] = df["Top_Genes"].map(clean_text).replace("", "Unassigned")
    df["Row_Label"] = df.apply(lambda r: f"{r['Gene_Label']} ({int(r['Max_Unique_Patients'])} pts)", axis=1)
    df["Tier_Num"] = df["Evidence_Tier"].map(tier_order)
    df["Bubble_Size"] = scale_sizes(df["Max_Unique_Patients"])

    save_tsv(df, outdir / "gwas_hotspot_besthit_dotplot_table.tsv")

    plot_df = df.iloc[::-1].copy()
    plot_df["y"] = range(len(plot_df))

    fig, ax = plt.subplots(figsize=(14, 9))
    sc = ax.scatter(
        plot_df["Distance_to_GWAS_bp"],
        plot_df["y"],
        s=plot_df["Bubble_Size"],
        c=plot_df["Tier_Num"],
        alpha=0.9,
    )

    ax.set_yticks(plot_df["y"])
    ax.set_yticklabels(plot_df["Row_Label"], fontsize=10)
    ax.set_xlabel("Distance to nearest CRC GWAS SNP (bp)")
    ax.set_ylabel("Hotspot linked gene")
    ax.set_title("CRC GWAS support for recurrent disrupted hotspots")
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Evidence tier rank")

    for _, row in plot_df.iterrows():
        ax.text(
            row["Distance_to_GWAS_bp"],
            row["y"] + 0.22,
            clean_text(row["GWAS_SNP"]),
            ha="center",
            fontsize=8,
        )

    fig.tight_layout()
    save_fig(fig, outdir / "gwas_hotspot_besthit_dotplot.png", dpi)


def gwas_overlap_schematic(
    gwas_besthit_file: Path | None,
    outdir: Path,
    dpi: int,
    top_n: int,
    require_gene_concordance: bool,
) -> None:
    df = maybe_read(gwas_besthit_file)

    if df.empty:
        print("[SKIP] GWAS schematic: no input file")
        return

    required = [
        "Hotspot_ID",
        "Chr",
        "Hotspot_Start",
        "Hotspot_End",
        "Hotspot_Width_bp",
        "Max_Unique_Patients",
        "Top_TFs",
        "Top_Genes",
        "GWAS_SNP",
        "GWAS_Pos",
        "Mapped_Gene",
        "Distance_to_GWAS_bp",
        "Evidence_Tier",
        "Gene_Concordance",
    ]
    if not check_columns(df, required, "GWAS schematic"):
        return

    df = df.copy()
    for col in ["Hotspot_Start", "Hotspot_End", "Hotspot_Width_bp", "Max_Unique_Patients", "GWAS_Pos", "Distance_to_GWAS_bp"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Hotspot_Start", "Hotspot_End", "GWAS_Pos", "Distance_to_GWAS_bp"]).copy()
    df = df[df["Evidence_Tier"].isin({"Tier 1: exact", "Tier 1: <=1kb"})].copy()

    if require_gene_concordance:
        df = df[df["Gene_Concordance"].astype(str).str.lower() == "yes"].copy()

    if "Trait_Priority" in df.columns:
        df["Trait_Priority"] = pd.to_numeric(df["Trait_Priority"], errors="coerce")
        df = df[df["Trait_Priority"].fillna(99) <= 1].copy()

    if df.empty:
        print("[SKIP] GWAS schematic: no rows after filters")
        return

    df = df.sort_values(["Distance_to_GWAS_bp", "Max_Unique_Patients"], ascending=[True, False])
    df = df.drop_duplicates(subset=["Hotspot_ID"]).head(top_n).copy()

    save_tsv(df, outdir / "crc_gwas_hotspot_direct_overlap_schematic_input.tsv")

    n = len(df)
    fig_height = max(6, n * 0.9 + 1.8)
    fig, ax = plt.subplots(figsize=(15, fig_height))

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.8, n - 0.2)
    ax.axis("off")

    ax.text(
        0.02,
        n - 0.05,
        "Recurrent hotspot",
        fontsize=12,
        fontweight="bold",
        va="bottom",
    )
    ax.text(
        0.72,
        n - 0.05,
        "CRC GWAS support",
        fontsize=12,
        fontweight="bold",
        va="bottom",
    )

    y_positions = np.arange(n)[::-1]

    for y, (_, row) in zip(y_positions, df.iterrows()):
        gene = split_items(row["Top_Genes"])[0] if split_items(row["Top_Genes"]) else clean_text(row["Mapped_Gene"])
        tf = ", ".join(split_items(row["Top_TFs"])[:3]) or "NA"

        ax.text(
            0.02,
            y,
            f"{gene}\n{row['Chr']}:{int(row['Hotspot_Start']):,}-{int(row['Hotspot_End']):,}\nTF: {tf}",
            fontsize=9,
            va="center",
        )

        local_min = min(int(row["Hotspot_Start"]), int(row["Hotspot_End"]), int(row["GWAS_Pos"]))
        local_max = max(int(row["Hotspot_Start"]), int(row["Hotspot_End"]), int(row["GWAS_Pos"]))
        span = max(local_max - local_min, 1)
        pad = max(20, min(150, int(span * 0.8)))
        local_min -= pad
        local_max += pad
        span = max(local_max - local_min, 1)

        x0 = 0.34
        x1 = 0.68

        def map_x(pos):
            return x0 + (pos - local_min) / span * (x1 - x0)

        hs_x0 = map_x(int(row["Hotspot_Start"]))
        hs_x1 = map_x(int(row["Hotspot_End"]))
        snp_x = map_x(int(row["GWAS_Pos"]))

        ax.plot([x0, x1], [y, y], color="black", linewidth=1)
        ax.add_patch(Rectangle((hs_x0, y - 0.12), max(hs_x1 - hs_x0, 0.005), 0.24, alpha=0.6))
        ax.plot([snp_x, snp_x], [y - 0.22, y + 0.22], color="black", linewidth=2)

        ax.text(
            0.72,
            y,
            f"{row['GWAS_SNP']}\n{row['Evidence_Tier']}\n{int(row['Distance_to_GWAS_bp'])} bp",
            fontsize=9,
            va="center",
        )

    fig.tight_layout()
    save_fig(fig, outdir / "crc_gwas_hotspot_direct_overlap_schematic.png", dpi)


def tf_prioritization_figure(disruptions_file: Path | None, outdir: Path, dpi: int, top_n: int) -> Path | None:
    df = maybe_read(disruptions_file)

    if df.empty:
        print("[SKIP] TF prioritization: no input file")
        return None

    required = ["Name", "Tumor_Sample_Barcode", "Chr", "Motif_Start", "Motif_End", "Entropy"]
    if not check_columns(df, required, "TF prioritization"):
        return None

    df = df.copy()
    df["TF"] = df["Name"].map(motif_to_tf)
    df = df[df["TF"] != ""].copy()

    locus_cols = ["Chr", "Motif_Start", "Motif_End"]
    variant_cols = [
        col for col in [
            "Variant_Chr",
            "Variant_Start",
            "Variant_End",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode",
        ]
        if col in df.columns
    ]

    tf_patients = df.groupby("TF")["Tumor_Sample_Barcode"].nunique().reset_index(name="Unique_Patients")
    tf_loci = df[["TF"] + locus_cols].drop_duplicates().groupby("TF").size().reset_index(name="Unique_Disrupted_Loci")
    tf_rows = df.groupby("TF").size().reset_index(name="Num_Total_Disrupting_Rows")
    tf_entropy = df.groupby("TF")["Entropy"].agg(["max", "mean"]).reset_index().rename(columns={"max": "Max_Entropy", "mean": "Mean_Entropy"})

    summary = tf_patients.merge(tf_loci, on="TF", how="outer")
    summary = summary.merge(tf_rows, on="TF", how="outer")
    summary = summary.merge(tf_entropy, on="TF", how="outer")

    if variant_cols:
        tf_variants = df[["TF"] + variant_cols].drop_duplicates().groupby("TF").size().reset_index(name="Unique_Variants")
        summary = summary.merge(tf_variants, on="TF", how="outer")
    else:
        summary["Unique_Variants"] = 0

    summary = summary.fillna(0)
    for col in ["Unique_Patients", "Unique_Disrupted_Loci", "Num_Total_Disrupting_Rows", "Unique_Variants"]:
        summary[col] = summary[col].astype(int)

    summary = summary.sort_values(
        ["Unique_Patients", "Unique_Disrupted_Loci", "Unique_Variants", "Max_Entropy"],
        ascending=[False, False, False, False],
    )

    out_table = outdir / "tf_prioritization.tsv"
    save_tsv(summary, out_table)

    plot_df = summary.head(top_n).sort_values("Unique_Patients", ascending=True)

    fig, ax = plt.subplots(figsize=(11, 8))
    ax.barh(plot_df["TF"], plot_df["Unique_Patients"])
    ax.set_xlabel("Unique patients")
    ax.set_ylabel("TF")
    ax.set_title("Top disrupted TF motifs by unique patients")

    maxv = plot_df["Unique_Patients"].max()
    for i, v in enumerate(plot_df["Unique_Patients"]):
        ax.text(v + maxv * 0.01, i, str(int(v)), va="center")

    fig.tight_layout()
    save_fig(fig, outdir / "top_tfs_by_unique_patients_detailed.png", dpi)

    return out_table


def tf_gene_intersection_bubble(intersection_file: Path | None, outdir: Path, dpi: int, top_tfs: int, top_genes: int) -> None:
    df = maybe_read(intersection_file)

    if df.empty:
        print("[SKIP] TF-gene bubble plot: no input file")
        return

    required = ["Chr", "Motif_Start", "Motif_End", "Unique_Patients_Locus", "Matching_Top_Genes", "Matching_Top_TFs"]
    if not check_columns(df, required, "TF-gene bubble plot"):
        return

    df["Locus_ID"] = (
        df["Chr"].astype(str)
        + ":"
        + df["Motif_Start"].astype(int).astype(str)
        + "-"
        + df["Motif_End"].astype(int).astype(str)
    )

    rows = []

    for _, row in df.iterrows():
        genes = split_items(row["Matching_Top_Genes"])
        tfs = split_items(row["Matching_Top_TFs"])
        for tf in tfs:
            for gene in genes:
                rows.append(
                    {
                        "TF": tf,
                        "Gene": gene,
                        "Locus_ID": row["Locus_ID"],
                        "Unique_Patients_Locus": int(row["Unique_Patients_Locus"]),
                    }
                )

    if not rows:
        print("[SKIP] TF-gene bubble plot: no TF-gene pairs")
        return

    pairs = pd.DataFrame(rows).drop_duplicates()

    summary = (
        pairs.groupby(["TF", "Gene"], as_index=False)
        .agg(
            Num_Recurrent_Loci=("Locus_ID", "nunique"),
            Max_Unique_Patients=("Unique_Patients_Locus", "max"),
            Mean_Unique_Patients=("Unique_Patients_Locus", "mean"),
        )
        .sort_values(["Max_Unique_Patients", "Num_Recurrent_Loci"], ascending=[False, False])
    )

    top_tf_list = (
        pairs.groupby("TF")["Unique_Patients_Locus"]
        .max()
        .sort_values(ascending=False)
        .head(top_tfs)
        .index.tolist()
    )

    top_gene_list = (
        pairs.groupby("Gene")["Unique_Patients_Locus"]
        .max()
        .sort_values(ascending=False)
        .head(top_genes)
        .index.tolist()
    )

    plot_df = summary[summary["TF"].isin(top_tf_list) & summary["Gene"].isin(top_gene_list)].copy()

    if plot_df.empty:
        print("[SKIP] TF-gene bubble plot: empty after top filtering")
        return

    save_tsv(plot_df, outdir / "tf_gene_intersection_bubble_plot_table.tsv")

    tf_order = (
        plot_df.groupby("TF")["Max_Unique_Patients"]
        .max()
        .sort_values(ascending=True)
        .index.tolist()
    )
    gene_order = (
        plot_df.groupby("Gene")["Max_Unique_Patients"]
        .max()
        .sort_values(ascending=True)
        .index.tolist()
    )

    tf_to_y = {tf: i for i, tf in enumerate(tf_order)}
    gene_to_x = {gene: i for i, gene in enumerate(gene_order)}

    plot_df["x"] = plot_df["Gene"].map(gene_to_x)
    plot_df["y"] = plot_df["TF"].map(tf_to_y)
    plot_df["Bubble_Size"] = scale_sizes(plot_df["Num_Recurrent_Loci"])

    fig, ax = plt.subplots(figsize=(max(10, 0.7 * len(gene_order) + 4), max(7, 0.55 * len(tf_order) + 3)))

    sc = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=plot_df["Bubble_Size"],
        c=plot_df["Max_Unique_Patients"],
        alpha=0.85,
    )

    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=45, ha="right")
    ax.set_yticks(range(len(tf_order)))
    ax.set_yticklabels(tf_order)
    ax.set_xlabel("Gene")
    ax.set_ylabel("TF")
    ax.set_title("TF-gene recurrent locus intersections")
    ax.grid(alpha=0.2)
    ax.set_axisbelow(True)

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Maximum unique patients")

    fig.tight_layout()
    save_fig(fig, outdir / "tf_gene_intersection_bubble_plot.png", dpi)


def pathway_dotplots(external_dir: Path, outdir: Path, dpi: int, top_n: int) -> None:
    for label, title in [
        ("top_genes", "Pathway enrichment for top recurrent genes"),
        ("ncg_genes", "Pathway enrichment for recurrent genes overlapping NCG"),
    ]:
        path = external_dir / f"{label}_pathway_table.tsv"
        df = maybe_read(path)

        if df.empty:
            print(f"[SKIP] pathway dotplot {label}: no table")
            continue

        required = ["Term", "Minus_log10_Adjusted_P", "Num_Genes"]
        if not check_columns(df, required, f"pathway dotplot {label}"):
            continue

        plot_df = df.head(top_n).iloc[::-1].copy()
        plot_df["BubbleSize"] = scale_sizes(plot_df["Num_Genes"])

        fig, ax = plt.subplots(figsize=(12, max(6, 0.55 * len(plot_df) + 2)))
        ax.scatter(
            plot_df["Minus_log10_Adjusted_P"],
            range(len(plot_df)),
            s=plot_df["BubbleSize"],
            alpha=0.85,
        )

        ax.set_yticks(range(len(plot_df)))
        ax.set_yticklabels(plot_df["Term"])
        ax.set_xlabel("-log10 adjusted p-value")
        ax.set_ylabel("Pathway")
        ax.set_title(title)
        ax.grid(axis="x", alpha=0.25)
        ax.set_axisbelow(True)

        for i, (_, row) in enumerate(plot_df.iterrows()):
            ax.text(
                row["Minus_log10_Adjusted_P"],
                i,
                str(int(row["Num_Genes"])),
                ha="center",
                va="center",
                fontsize=8,
                fontweight="bold",
            )

        fig.tight_layout()
        save_fig(fig, outdir / f"{label}_pathway_dotplot.png", dpi)


def run_scoring(args) -> None:
    scoring_summary_figure(args.total_motifs, args.functional_motifs, args.outdir, args.dpi)

    tf_background = args.tf_background
    if tf_background is None and args.functional_table and args.disruptions:
        tf_background = tf_background_comparison(args.functional_table, args.disruptions, args.outdir)

    top_normalized_tfs_figure(tf_background, args.outdir, args.dpi, args.top_tfs)


def run_hotspots(args) -> None:
    recurrent_loci_manhattan(
        args.recurrent_dir / "recurrent_loci_full.tsv",
        args.outdir,
        args.dpi,
        args.min_patients_plot,
        args.top_labels,
    )

    hotspot_landscape_figure(
        args.recurrent_dir / "recurrent_hotspot_regions.tsv",
        args.outdir,
        args.dpi,
        args.top_labels,
    )


def run_loci(args) -> None:
    top_integrated_loci_figure(
        args.recurrent_dir / "integrated_recurrent_locus_table.tsv",
        args.outdir,
        args.dpi,
        args.top_loci,
    )


def run_gwas(args) -> None:
    besthit = args.external_dir / "crc_gwas_hotspot_besthit_shortlist.tsv"

    if not besthit.exists():
        besthit = args.external_dir / "crc_gwas_hotspot_besthit.tsv"

    gwas_hotspot_dotplot(besthit, args.outdir, args.dpi, args.top_gwas)
    gwas_overlap_schematic(
        besthit,
        args.outdir,
        args.dpi,
        args.top_gwas_schematic,
        args.require_gene_concordance,
    )


def run_tf_gene(args) -> None:
    tf_prioritization_figure(args.disruptions, args.outdir, args.dpi, args.top_tfs)

    tf_gene_intersection_bubble(
        args.recurrent_dir / "top_gene_tf_locus_intersections.tsv",
        args.outdir,
        args.dpi,
        args.top_tfs,
        args.top_genes,
    )


def run_pathways(args) -> None:
    pathway_dotplots(args.external_dir, args.outdir, args.dpi, args.top_pathways)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate thesis-ready figures from final analysis tables."
    )

    parser.add_argument(
        "--group",
        choices=["all", "scoring", "hotspots", "loci", "gwas", "tf_gene", "pathways"],
        default="all",
        help="Figure group to generate."
    )
    parser.add_argument(
        "--recurrent-dir",
        type=Path,
        default=Path("results/tables/recurrent_hotspots_and_genes"),
        help="Directory from 06_recurrent_hotspots_and_genes.py."
    )
    parser.add_argument(
        "--external-dir",
        type=Path,
        default=Path("results/tables/external_interpretation"),
        help="Directory from 07_external_interpretation.py."
    )
    parser.add_argument(
        "--variant-dir",
        type=Path,
        default=Path("results/tables/variant_overlap_disruption"),
        help="Directory from 05_variant_overlap_and_disruption.py."
    )
    parser.add_argument(
        "--disruptions",
        type=Path,
        default=None,
        help="Motif disrupting variants table."
    )
    parser.add_argument(
        "--functional-table",
        type=Path,
        default=None,
        help="Functional motif prediction table, used for normalized TF background."
    )
    parser.add_argument(
        "--tf-background",
        type=Path,
        default=None,
        help="Optional precomputed TF background comparison table."
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/figures"),
        help="Output figure directory."
    )
    parser.add_argument("--dpi", type=int, default=300)

    parser.add_argument("--total-motifs", type=int, default=0)
    parser.add_argument("--functional-motifs", type=int, default=0)

    parser.add_argument("--top-loci", type=int, default=15)
    parser.add_argument("--top-labels", type=int, default=8)
    parser.add_argument("--top-gwas", type=int, default=15)
    parser.add_argument("--top-gwas-schematic", type=int, default=10)
    parser.add_argument("--top-tfs", type=int, default=20)
    parser.add_argument("--top-genes", type=int, default=12)
    parser.add_argument("--top-pathways", type=int, default=12)
    parser.add_argument("--min-patients-plot", type=int, default=2)
    parser.add_argument("--require-gene-concordance", action="store_true")

    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.disruptions is None:
        candidate = args.variant_dir / "motif_disrupting_variants.tsv"
        if candidate.exists():
            args.disruptions = candidate

    selected = (
        ["scoring", "hotspots", "loci", "gwas", "tf_gene", "pathways"]
        if args.group == "all"
        else [args.group]
    )

    if "scoring" in selected:
        run_scoring(args)

    if "hotspots" in selected:
        run_hotspots(args)

    if "loci" in selected:
        run_loci(args)

    if "gwas" in selected:
        run_gwas(args)

    if "tf_gene" in selected:
        run_tf_gene(args)

    if "pathways" in selected:
        run_pathways(args)

    print(f"Done. Figures saved in: {args.outdir}")


if __name__ == "__main__":
    main()
