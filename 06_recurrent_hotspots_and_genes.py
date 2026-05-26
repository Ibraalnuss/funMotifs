#!/usr/bin/env python3

"""
06_recurrent_hotspots_and_genes.py

Identify recurrent disrupted loci and build gene/TF prioritization tables
from the motif-disrupting variant output of script 05.

Steps:
  1. Count unique patients per disrupted locus; filter by --min-patients.
  2. Merge nearby recurrent loci into hotspot regions (--merge-gap).
  3. Build gene and TF prioritization tables linked to recurrent loci.
  4. Produce an integrated locus/TF/gene summary table.
  5. Optionally annotate hotspots with functional motif feature summaries
     from the ROSETTA input chunks.
  6. Optionally extract non-TSS hotspot examples for thesis interpretation.

Note: GWAS and NCG integration is handled in 07_external_interpretation.py,
not here. This script focuses on recurrence patterns within the CRC cohort.

Main outputs:
    recurrent_loci_full.tsv
    recurrent_loci_filtered.tsv
    integrated_recurrent_locus_table.tsv
    gene_prioritization_from_recurrent_loci.tsv
    recurrent_hotspot_regions.tsv
    recurrent_hotspot_regions.bed
    gene_to_top_tfs.tsv
    tf_to_top_genes.tsv
    top_hotspot_functional_motif_summary.tsv      (optional)
    non_tss_hotspot_examples.tsv                  (optional)
    recurrent_hotspots_and_genes_summary.txt

Example:
    python scripts/06_recurrent_hotspots_and_genes.py \\
      --disruptions results/tables/variant_overlap_disruption/motif_disrupting_variants.tsv \\
      --outdir results/tables/recurrent_hotspots_and_genes \\
      --min-patients 10 \\
      --merge-gap 200
"""

from __future__ import annotations

import argparse
import glob
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


LOCUS_COLS = ["Chr", "Motif_Start", "Motif_End"]
VARIANT_COLS = [
    "Variant_Chr",
    "Variant_Start",
    "Variant_End",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
]

DEFAULT_TOP_GENES_OF_INTEREST = {
    "SRPK1", "CUX1", "MAPK1", "PRKDC", "FOXP1", "PREX1", "CDK12", "TRIO"
}

DEFAULT_TOP_TFS_OF_INTEREST = {
    "CTCF", "KLF12", "KLF1", "KLF2", "KLF3", "KLF9", "KLF11", "KLF14",
    "FOXD3", "RREB1", "ZNF384", "ZNF460"
}

FEATURE_COLS = [
    "CTCF.only_CTCF.bound", "DNase.H3K4me3", "DNase.H3K4me3_CTCF.bound",
    "DNase.only", "Low.DNase", "PLS", "PLS_CTCF.bound",
    "dELS", "dELS_CTCF.bound", "pELS", "pELS_CTCF.bound",
    "Enh", "Quies", "Repr", "TSS", "dnase__seq", "fantom",
    "footprints", "DTZ", "ERD", "LRD", "UTZ",
    "numothertfbinding_wang", "tfexpr_wang"
]


def clean_text(value) -> str:
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value.lower() in {"", "nan", "na", "none"}:
        return ""
    return " ".join(value.split())


def motif_to_tf(motif_name: str) -> str:
    motif_name = clean_text(motif_name)
    if not motif_name:
        return ""
    if "_" in motif_name:
        return motif_name.rsplit("_", 1)[0].strip()
    return motif_name


def split_items(text) -> list[str]:
    text = clean_text(text)
    if not text:
        return []
    return [x.strip() for x in text.split(",") if x.strip()]


def top_join(series, n: int) -> str:
    values = []
    seen = set()

    for item in series:
        for x in split_items(item):
            if x not in seen:
                seen.add(x)
                values.append(x)

    return ", ".join(values[:n])


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def save_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def load_gene_map(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}

    df = pd.read_csv(path, sep="\t", low_memory=False)

    gene_col = None
    symbol_col = None

    for candidate in ["Gene", "gene_id", "ensembl_gene_id"]:
        if candidate in df.columns:
            gene_col = candidate
            break

    for candidate in ["Gene_Symbol", "gene_symbol", "hgnc_symbol", "external_gene_name"]:
        if candidate in df.columns:
            symbol_col = candidate
            break

    if gene_col is None or symbol_col is None:
        raise ValueError(
            f"Gene map must contain a gene column and a symbol column: {path}"
        )

    out = {}
    for _, row in df[[gene_col, symbol_col]].iterrows():
        gene = clean_text(row[gene_col])
        symbol = clean_text(row[symbol_col])
        if gene and symbol:
            out[gene] = symbol

    return out


def map_gene_label(gene_value, gene_map: dict[str, str]) -> str:
    parts = split_items(gene_value)
    labels = []

    for gene in parts:
        mapped = gene_map.get(gene, gene)
        mapped = clean_text(mapped)
        if mapped:
            labels.append(mapped)

    seen = set()
    out = []

    for label in labels:
        if label not in seen:
            seen.add(label)
            out.append(label)

    return ",".join(out)


def add_tf_and_gene_labels(df: pd.DataFrame, gene_map: dict[str, str]) -> pd.DataFrame:
    df = df.copy()
    df["TF"] = df["Name"].map(motif_to_tf)
    df["Gene_Label"] = df["Gene"].apply(lambda x: map_gene_label(x, gene_map))
    return df


def build_recurrent_loci(df: pd.DataFrame, min_patients: int, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    locus_patient = df[LOCUS_COLS + ["Tumor_Sample_Barcode"]].drop_duplicates()

    recurrent_loci = (
        locus_patient
        .groupby(LOCUS_COLS)["Tumor_Sample_Barcode"]
        .nunique()
        .reset_index(name="Unique_Patients")
        .sort_values("Unique_Patients", ascending=False)
        .reset_index(drop=True)
    )

    recurrent_loci_filtered = recurrent_loci[
        recurrent_loci["Unique_Patients"] >= min_patients
    ].copy()

    recurrent_loci_filtered = recurrent_loci_filtered.rename(
        columns={"Unique_Patients": "Locus_Unique_Patients"}
    )

    save_tsv(recurrent_loci, outdir / "recurrent_loci_full.tsv")
    save_tsv(recurrent_loci_filtered, outdir / "recurrent_loci_filtered.tsv")

    return recurrent_loci, recurrent_loci_filtered


def build_gene_prioritization(
    df: pd.DataFrame,
    recurrent_loci_filtered: pd.DataFrame,
    outdir: Path,
) -> pd.DataFrame:
    if recurrent_loci_filtered.empty:
        gene_summary = pd.DataFrame()
        save_tsv(gene_summary, outdir / "gene_prioritization_from_recurrent_loci.tsv")
        return gene_summary

    merged = df.merge(
        recurrent_loci_filtered,
        on=LOCUS_COLS,
        how="inner",
    )

    merged = merged[merged["Gene"].map(clean_text) != ""].copy()

    if merged.empty:
        gene_summary = pd.DataFrame()
        save_tsv(gene_summary, outdir / "gene_prioritization_from_recurrent_loci.tsv")
        return gene_summary

    gene_patients = (
        merged.groupby("Gene")["Tumor_Sample_Barcode"]
        .nunique()
        .reset_index(name="Unique_Patients")
    )

    gene_loci = (
        merged[["Gene"] + LOCUS_COLS]
        .drop_duplicates()
        .groupby("Gene")
        .size()
        .reset_index(name="Unique_Recurrent_Loci")
    )

    gene_variants = (
        merged[["Gene"] + VARIANT_COLS]
        .drop_duplicates()
        .groupby("Gene")
        .size()
        .reset_index(name="Unique_Variants")
    )

    gene_rows = (
        merged.groupby("Gene")
        .size()
        .reset_index(name="Num_Total_Disrupting_Rows")
    )

    gene_max_locus_patients = (
        merged.groupby("Gene")["Locus_Unique_Patients"]
        .max()
        .reset_index(name="Max_Locus_Unique_Patients")
    )

    gene_entropy = (
        merged.groupby("Gene")["Entropy"]
        .agg(["max", "mean"])
        .reset_index()
        .rename(columns={"max": "Max_Entropy", "mean": "Mean_Entropy"})
    )

    gene_summary = gene_patients.merge(gene_loci, on="Gene", how="outer")
    gene_summary = gene_summary.merge(gene_variants, on="Gene", how="outer")
    gene_summary = gene_summary.merge(gene_rows, on="Gene", how="outer")
    gene_summary = gene_summary.merge(gene_max_locus_patients, on="Gene", how="outer")
    gene_summary = gene_summary.merge(gene_entropy, on="Gene", how="outer")
    gene_summary = gene_summary.fillna(0)

    int_cols = [
        "Unique_Patients",
        "Unique_Recurrent_Loci",
        "Unique_Variants",
        "Num_Total_Disrupting_Rows",
        "Max_Locus_Unique_Patients",
    ]

    for col in int_cols:
        gene_summary[col] = gene_summary[col].astype(int)

    label_map = (
        merged[["Gene", "Gene_Label"]]
        .drop_duplicates()
        .set_index("Gene")["Gene_Label"]
        .to_dict()
    )

    gene_summary["Gene_Label"] = gene_summary["Gene"].map(label_map).fillna(gene_summary["Gene"])

    gene_summary = gene_summary[
        [
            "Gene_Label",
            "Gene",
            "Unique_Patients",
            "Unique_Recurrent_Loci",
            "Unique_Variants",
            "Num_Total_Disrupting_Rows",
            "Max_Locus_Unique_Patients",
            "Max_Entropy",
            "Mean_Entropy",
        ]
    ]

    gene_summary = gene_summary.sort_values(
        [
            "Unique_Patients",
            "Unique_Recurrent_Loci",
            "Max_Locus_Unique_Patients",
            "Unique_Variants",
            "Max_Entropy",
        ],
        ascending=[False, False, False, False, False],
    ).reset_index(drop=True)

    save_tsv(gene_summary, outdir / "gene_prioritization_from_recurrent_loci.tsv")
    save_tsv(
        gene_summary[
            [
                "Gene_Label",
                "Gene",
                "Unique_Patients",
                "Unique_Recurrent_Loci",
                "Max_Locus_Unique_Patients",
                "Unique_Variants",
                "Max_Entropy",
            ]
        ],
        outdir / "pathway_input_genes_from_recurrent_loci.tsv",
    )

    return gene_summary


def build_integrated_locus_tables(
    df: pd.DataFrame,
    recurrent_loci_filtered: pd.DataFrame,
    outdir: Path,
    top_n: int,
    genes_of_interest: set[str],
    tfs_of_interest: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if recurrent_loci_filtered.empty:
        empty = pd.DataFrame()
        save_tsv(empty, outdir / "integrated_recurrent_locus_table.tsv")
        save_tsv(empty, outdir / "gene_to_top_tfs.tsv")
        save_tsv(empty, outdir / "tf_to_top_genes.tsv")
        save_tsv(empty, outdir / "top_gene_tf_locus_intersections.tsv")
        return empty, empty, empty, empty

    locus_summary = (
        df.groupby(LOCUS_COLS, as_index=False)
        .agg(
            Unique_Patients_Locus=("Tumor_Sample_Barcode", "nunique"),
            Num_Disrupting_Rows_Locus=("Tumor_Sample_Barcode", "size"),
            Max_Entropy_Locus=("Entropy", "max"),
            Mean_Entropy_Locus=("Entropy", "mean"),
        )
    )

    unique_variants = (
        df[LOCUS_COLS + VARIANT_COLS]
        .drop_duplicates()
        .groupby(LOCUS_COLS, as_index=False)
        .size()
        .rename(columns={"size": "Unique_Variants_Locus"})
    )

    locus_summary = locus_summary.merge(unique_variants, on=LOCUS_COLS, how="left")
    locus_summary["Unique_Variants_Locus"] = (
        locus_summary["Unique_Variants_Locus"].fillna(0).astype(int)
    )

    recurrent = locus_summary.merge(
        recurrent_loci_filtered[LOCUS_COLS],
        on=LOCUS_COLS,
        how="inner",
    )

    recurrent_rows = df.merge(recurrent[LOCUS_COLS], on=LOCUS_COLS, how="inner")

    tf_locus = (
        recurrent_rows[recurrent_rows["TF"] != ""]
        .groupby(LOCUS_COLS + ["TF"], as_index=False)
        .agg(Unique_Patients_TF=("Tumor_Sample_Barcode", "nunique"))
    )

    if not tf_locus.empty:
        tf_per_locus = (
            tf_locus.groupby(LOCUS_COLS, as_index=False)
            .apply(
                lambda g: pd.Series(
                    {
                        "Top_TFs": ", ".join(
                            g.sort_values(
                                ["Unique_Patients_TF", "TF"],
                                ascending=[False, True],
                            )["TF"].head(top_n).astype(str)
                        ),
                        "Num_Unique_TFs": g["TF"].nunique(),
                    }
                )
            )
            .reset_index(drop=True)
        )
    else:
        tf_per_locus = pd.DataFrame(columns=LOCUS_COLS + ["Top_TFs", "Num_Unique_TFs"])

    gene_long_rows = []

    for _, row in recurrent_rows.iterrows():
        for gene in split_items(row["Gene_Label"]):
            gene_long_rows.append(
                {
                    "Chr": row["Chr"],
                    "Motif_Start": row["Motif_Start"],
                    "Motif_End": row["Motif_End"],
                    "Gene": gene,
                    "Tumor_Sample_Barcode": row["Tumor_Sample_Barcode"],
                }
            )

    gene_long = pd.DataFrame(gene_long_rows)

    if not gene_long.empty:
        gene_long = gene_long.drop_duplicates()

        gene_locus = (
            gene_long.groupby(LOCUS_COLS + ["Gene"], as_index=False)
            .agg(Unique_Patients_Gene=("Tumor_Sample_Barcode", "nunique"))
        )

        gene_per_locus = (
            gene_locus.groupby(LOCUS_COLS, as_index=False)
            .apply(
                lambda g: pd.Series(
                    {
                        "Top_Genes": ", ".join(
                            g.sort_values(
                                ["Unique_Patients_Gene", "Gene"],
                                ascending=[False, True],
                            )["Gene"].head(top_n).astype(str)
                        ),
                        "Num_Unique_Genes": g["Gene"].nunique(),
                    }
                )
            )
            .reset_index(drop=True)
        )
    else:
        gene_per_locus = pd.DataFrame(columns=LOCUS_COLS + ["Top_Genes", "Num_Unique_Genes"])

    integrated = recurrent.merge(tf_per_locus, on=LOCUS_COLS, how="left")
    integrated = integrated.merge(gene_per_locus, on=LOCUS_COLS, how="left")

    integrated["Top_TFs"] = integrated["Top_TFs"].fillna("")
    integrated["Top_Genes"] = integrated["Top_Genes"].fillna("")
    integrated["Num_Unique_TFs"] = integrated["Num_Unique_TFs"].fillna(0).astype(int)
    integrated["Num_Unique_Genes"] = integrated["Num_Unique_Genes"].fillna(0).astype(int)

    integrated = integrated.sort_values(
        ["Unique_Patients_Locus", "Unique_Variants_Locus", "Num_Disrupting_Rows_Locus"],
        ascending=[False, False, False],
    ).reset_index(drop=True)

    save_tsv(integrated, outdir / "integrated_recurrent_locus_table.tsv")

    gene_tf_table, tf_gene_table = build_gene_tf_tables(recurrent_rows, gene_long, outdir, top_n)
    intersection = build_interest_intersection_table(
        integrated,
        genes_of_interest,
        tfs_of_interest,
        outdir,
    )

    return integrated, gene_tf_table, tf_gene_table, intersection


def build_gene_tf_tables(
    recurrent_rows: pd.DataFrame,
    gene_long: pd.DataFrame,
    outdir: Path,
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if gene_long.empty:
        gene_tf_table = pd.DataFrame(columns=["Gene", "Top_TFs", "Max_TF_Unique_Patients", "Num_TFs_Linked"])
        tf_gene_table = pd.DataFrame(columns=["TF", "Top_Genes", "Max_Gene_Unique_Patients", "Num_Genes_Linked"])
        save_tsv(gene_tf_table, outdir / "gene_to_top_tfs.tsv")
        save_tsv(tf_gene_table, outdir / "tf_to_top_genes.tsv")
        return gene_tf_table, tf_gene_table

    gene_tf = (
        recurrent_rows[recurrent_rows["TF"] != ""]
        .merge(
            gene_long.drop_duplicates(),
            on=["Chr", "Motif_Start", "Motif_End", "Tumor_Sample_Barcode"],
            how="inner",
        )
    )

    gene_tf_summary = (
        gene_tf.groupby(["Gene", "TF"], as_index=False)
        .agg(
            Unique_Patients=("Tumor_Sample_Barcode", "nunique"),
            Unique_Loci=("Chr", "count"),
        )
    )

    gene_tf_table = (
        gene_tf_summary.groupby("Gene", as_index=False)
        .apply(
            lambda g: pd.Series(
                {
                    "Top_TFs": ", ".join(
                        g.sort_values(["Unique_Patients", "TF"], ascending=[False, True])["TF"]
                        .head(top_n)
                        .astype(str)
                    ),
                    "Max_TF_Unique_Patients": int(g["Unique_Patients"].max()),
                    "Num_TFs_Linked": g["TF"].nunique(),
                }
            )
        )
        .reset_index(drop=True)
        .sort_values("Max_TF_Unique_Patients", ascending=False)
        .reset_index(drop=True)
    )

    tf_gene_table = (
        gene_tf_summary.groupby("TF", as_index=False)
        .apply(
            lambda g: pd.Series(
                {
                    "Top_Genes": ", ".join(
                        g.sort_values(["Unique_Patients", "Gene"], ascending=[False, True])["Gene"]
                        .head(top_n)
                        .astype(str)
                    ),
                    "Max_Gene_Unique_Patients": int(g["Unique_Patients"].max()),
                    "Num_Genes_Linked": g["Gene"].nunique(),
                }
            )
        )
        .reset_index(drop=True)
        .sort_values("Max_Gene_Unique_Patients", ascending=False)
        .reset_index(drop=True)
    )

    save_tsv(gene_tf_table, outdir / "gene_to_top_tfs.tsv")
    save_tsv(tf_gene_table, outdir / "tf_to_top_genes.tsv")

    return gene_tf_table, tf_gene_table


def build_interest_intersection_table(
    integrated: pd.DataFrame,
    genes_of_interest: set[str],
    tfs_of_interest: set[str],
    outdir: Path,
) -> pd.DataFrame:
    rows = []

    for _, row in integrated.iterrows():
        genes = {x.strip() for x in str(row["Top_Genes"]).split(",") if x.strip()}
        tfs = {x.strip() for x in str(row["Top_TFs"]).split(",") if x.strip()}

        matching_genes = sorted(genes & genes_of_interest)
        matching_tfs = sorted(tfs & tfs_of_interest)

        if matching_genes and matching_tfs:
            rows.append(
                {
                    "Chr": row["Chr"],
                    "Motif_Start": row["Motif_Start"],
                    "Motif_End": row["Motif_End"],
                    "Unique_Patients_Locus": row["Unique_Patients_Locus"],
                    "Unique_Variants_Locus": row["Unique_Variants_Locus"],
                    "Num_Disrupting_Rows_Locus": row["Num_Disrupting_Rows_Locus"],
                    "Max_Entropy_Locus": row["Max_Entropy_Locus"],
                    "Matching_Top_Genes": ", ".join(matching_genes),
                    "Matching_Top_TFs": ", ".join(matching_tfs),
                    "Top_Genes": row["Top_Genes"],
                    "Top_TFs": row["Top_TFs"],
                }
            )

    out = pd.DataFrame(rows)

    if not out.empty:
        out = out.sort_values(
            ["Unique_Patients_Locus", "Unique_Variants_Locus"],
            ascending=[False, False],
        ).reset_index(drop=True)

    save_tsv(out, outdir / "top_gene_tf_locus_intersections.tsv")
    return out


def build_hotspots(
    integrated: pd.DataFrame,
    outdir: Path,
    merge_gap: int,
    top_n_tfs: int,
    top_n_genes: int,
) -> pd.DataFrame:
    if integrated.empty:
        out = pd.DataFrame()
        save_tsv(out, outdir / "recurrent_hotspot_regions.tsv")
        save_tsv(out, outdir / "recurrent_hotspot_regions.bed")
        return out

    df = integrated.copy()
    df["Chr"] = df["Chr"].map(clean_text)
    df["Motif_Start"] = df["Motif_Start"].astype(int)
    df["Motif_End"] = df["Motif_End"].astype(int)
    df = df.sort_values(["Chr", "Motif_Start", "Motif_End"]).reset_index(drop=True)

    hotspot_rows = []
    hotspot_id = 0

    for chr_value, sub in df.groupby("Chr", sort=False):
        sub = sub.sort_values(["Motif_Start", "Motif_End"]).reset_index(drop=True)

        current = []
        current_start = None
        current_end = None

        for _, row in sub.iterrows():
            start = int(row["Motif_Start"])
            end = int(row["Motif_End"])

            if not current:
                current = [row]
                current_start = start
                current_end = end
                continue

            if start <= current_end + merge_gap:
                current.append(row)
                current_end = max(current_end, end)
            else:
                hotspot_id += 1
                hotspot_rows.append(
                    summarize_hotspot(
                        hotspot_id,
                        chr_value,
                        current,
                        current_start,
                        current_end,
                        top_n_tfs,
                        top_n_genes,
                    )
                )
                current = [row]
                current_start = start
                current_end = end

        if current:
            hotspot_id += 1
            hotspot_rows.append(
                summarize_hotspot(
                    hotspot_id,
                    chr_value,
                    current,
                    current_start,
                    current_end,
                    top_n_tfs,
                    top_n_genes,
                )
            )

    out = pd.DataFrame(hotspot_rows)

    out = out.sort_values(
        ["Max_Unique_Patients", "Num_Loci_In_Hotspot", "Hotspot_Width_bp"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    save_tsv(out, outdir / "recurrent_hotspot_regions.tsv")

    bed = out[["Chr", "Hotspot_Start", "Hotspot_End", "Hotspot_ID"]].copy()
    save_tsv(bed, outdir / "recurrent_hotspot_regions.bed")

    return out


def summarize_hotspot(
    hotspot_id: int,
    chr_value: str,
    rows: list[pd.Series],
    start: int,
    end: int,
    top_n_tfs: int,
    top_n_genes: int,
) -> dict:
    df = pd.DataFrame(rows)

    representative = (
        df.sort_values(
            ["Unique_Patients_Locus", "Motif_Start", "Motif_End"],
            ascending=[False, True, True],
        )
        .apply(lambda r: f"{r['Chr']}:{int(r['Motif_Start'])}-{int(r['Motif_End'])}", axis=1)
        .iloc[0]
    )

    return {
        "Hotspot_ID": f"HS_{hotspot_id}",
        "Chr": chr_value,
        "Hotspot_Start": start,
        "Hotspot_End": end,
        "Hotspot_Width_bp": end - start + 1,
        "Num_Loci_In_Hotspot": len(df),
        "Max_Unique_Patients": int(df["Unique_Patients_Locus"].max()),
        "Mean_Unique_Patients": float(df["Unique_Patients_Locus"].mean()),
        "Top_TFs": top_join(df["Top_TFs"], top_n_tfs) if "Top_TFs" in df.columns else "",
        "Top_Genes": top_join(df["Top_Genes"], top_n_genes) if "Top_Genes" in df.columns else "",
        "Representative_Locus": representative,
    }


def summarize_event_string(group: pd.DataFrame) -> str:
    events = []

    grouped = (
        group.groupby(["Reference_Allele", "Tumor_Seq_Allele2"], dropna=False)["Tumor_Sample_Barcode"]
        .nunique()
        .reset_index(name="Unique_Samples")
        .sort_values("Unique_Samples", ascending=False)
    )

    for _, row in grouped.iterrows():
        events.append(
            f"{row['Reference_Allele']}->{row['Tumor_Seq_Allele2']} "
            f"({int(row['Unique_Samples'])} samples)"
        )

    return "; ".join(events)


def summarize_features(row: pd.Series) -> str:
    values = []

    for col in FEATURE_COLS:
        if col not in row.index:
            continue

        try:
            value = float(row[col])
        except Exception:
            continue

        if value != 0:
            values.append(f"{col}={value:g}")

    return "; ".join(values) if values else "all_zero"


def build_hotspot_motif_summary(
    hotspots: pd.DataFrame,
    disruption_df: pd.DataFrame,
    functional_bed: Path,
    functional_table: Path | None,
    chunk_glob: str | None,
    outdir: Path,
    top_n_hotspots: int,
) -> pd.DataFrame:
    if hotspots.empty:
        out = pd.DataFrame()
        save_tsv(out, outdir / "top_hotspot_functional_motif_summary.tsv")
        return out

    if functional_bed is None or not functional_bed.exists():
        print("[SKIP] No functional BED provided for hotspot motif summary.")
        return pd.DataFrame()

    selected_hotspots = hotspots.head(top_n_hotspots).copy()

    bed = pd.read_csv(
        functional_bed,
        sep="\t",
        header=None,
        names=["Motif_Chr", "Motif_Start", "Motif_End", "Motif_Name", "Strand", "mid"],
    )

    records = []

    for _, hotspot in selected_hotspots.iterrows():
        sub = disruption_df[
            (disruption_df["Chr"] == hotspot["Chr"])
            & (disruption_df["Motif_Start"] <= int(hotspot["Hotspot_End"]))
            & (disruption_df["Motif_End"] >= int(hotspot["Hotspot_Start"]))
        ].copy()

        if sub.empty:
            continue

        for keys, group in sub.groupby(["Chr", "Motif_Start", "Motif_End", "Name", "Strand"], dropna=False):
            records.append(
                {
                    "Hotspot_ID": hotspot["Hotspot_ID"],
                    "Top_Genes": hotspot.get("Top_Genes", ""),
                    "Top_TFs": hotspot.get("Top_TFs", ""),
                    "Representative_Locus": hotspot.get("Representative_Locus", ""),
                    "Motif_Chr": keys[0],
                    "Motif_Start": keys[1],
                    "Motif_End": keys[2],
                    "Motif_Name": keys[3],
                    "Strand": keys[4],
                    "Disruption_Events": summarize_event_string(group),
                }
            )

    motif_summary = pd.DataFrame(records).drop_duplicates()

    if motif_summary.empty:
        save_tsv(motif_summary, outdir / "top_hotspot_functional_motif_summary.tsv")
        return motif_summary

    merged = motif_summary.merge(
        bed,
        on=["Motif_Chr", "Motif_Start", "Motif_End", "Motif_Name", "Strand"],
        how="left",
    )

    if functional_table and functional_table.exists():
        fun = pd.read_csv(functional_table, sep="\t", low_memory=False)
        keep_cols = [c for c in ["mid", "predictedClass", "1", "0"] if c in fun.columns]
        if "mid" in keep_cols:
            merged = merged.merge(fun[keep_cols], on="mid", how="left")

    target_mids = set(merged["mid"].dropna().astype(int).tolist())

    if chunk_glob:
        feature_rows = []
        for path in sorted(glob.glob(chunk_glob)):
            chunk = pd.read_csv(path, sep="\t", low_memory=False)
            sub = chunk[chunk["mid"].isin(target_mids)].copy()
            if not sub.empty:
                sub["source_chunk"] = Path(path).name
                feature_rows.append(sub)

        if feature_rows:
            features = pd.concat(feature_rows, ignore_index=True)
            features["Functional_Feature_Summary"] = features.apply(summarize_features, axis=1)
            merged = merged.merge(
                features[["mid", "Functional_Feature_Summary"]].drop_duplicates(),
                on="mid",
                how="left",
            )

    sort_cols = [c for c in ["Hotspot_ID", "Motif_Start", "Motif_Name", "mid"] if c in merged.columns]
    if sort_cols:
        merged = merged.sort_values(sort_cols).drop_duplicates()

    save_tsv(merged, outdir / "top_hotspot_functional_motif_summary.tsv")
    return merged


def total_samples_from_events(text: str) -> int:
    nums = re.findall(r"\((\d+)\s+samples\)", str(text))
    return sum(int(x) for x in nums)


def build_non_tss_examples(
    motif_summary: pd.DataFrame,
    outdir: Path,
) -> pd.DataFrame:
    if motif_summary.empty or "Functional_Feature_Summary" not in motif_summary.columns:
        out = pd.DataFrame()
        save_tsv(out, outdir / "non_tss_hotspot_examples.tsv")
        return out

    support_terms = ["Low.DNase", "dnase__seq", "Enh", "dELS", "pELS", "footprints", "fantom"]
    pattern = "|".join(support_terms)

    non_tss = motif_summary[
        ~motif_summary["Functional_Feature_Summary"].fillna("").str.contains("TSS", case=False, na=False)
    ].copy()

    non_tss = non_tss[
        non_tss["Functional_Feature_Summary"].fillna("").str.contains(pattern, case=False, na=False)
    ].copy()

    if "Disruption_Events" in non_tss.columns:
        non_tss["Total_Disrupted_Samples"] = non_tss["Disruption_Events"].apply(total_samples_from_events)

    keep_cols = [
        "Hotspot_ID", "Top_Genes", "Top_TFs", "Representative_Locus",
        "Motif_Name", "mid", "Disruption_Events", "Functional_Feature_Summary",
        "Total_Disrupted_Samples",
    ]

    keep_cols = [c for c in keep_cols if c in non_tss.columns]
    non_tss = non_tss[keep_cols].sort_values(
        [c for c in ["Total_Disrupted_Samples", "Motif_Name"] if c in keep_cols],
        ascending=[False, True] if "Total_Disrupted_Samples" in keep_cols else True,
    )

    save_tsv(non_tss, outdir / "non_tss_hotspot_examples.tsv")
    return non_tss


def write_summary(
    outdir: Path,
    df: pd.DataFrame,
    recurrent_loci: pd.DataFrame,
    recurrent_loci_filtered: pd.DataFrame,
    gene_summary: pd.DataFrame,
    integrated: pd.DataFrame,
    hotspots: pd.DataFrame,
    motif_summary: pd.DataFrame,
    non_tss: pd.DataFrame,
    min_patients: int,
    merge_gap: int,
) -> None:
    lines = [
        "Recurrent hotspots and genes summary",
        "",
        f"Input disruption rows: {len(df):,}",
        f"All recurrent loci: {len(recurrent_loci):,}",
        f"Recurrent loci retained with >= {min_patients} patients: {len(recurrent_loci_filtered):,}",
        f"Genes prioritized: {len(gene_summary):,}",
        f"Integrated recurrent locus rows: {len(integrated):,}",
        f"Hotspots created: {len(hotspots):,}",
        f"Hotspot merge gap: {merge_gap} bp",
        f"Hotspot motif summary rows: {len(motif_summary):,}",
        f"Non-TSS hotspot examples: {len(non_tss):,}",
    ]

    if not hotspots.empty:
        lines.extend(
            [
                "",
                "Top 10 hotspots:",
                hotspots.head(10).to_string(index=False),
            ]
        )

    save_text(outdir / "recurrent_hotspots_and_genes_summary.txt", "\n".join(lines))


def run(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)

    if not args.disruptions.exists():
        raise FileNotFoundError(f"Disruption file not found: {args.disruptions}")

    print(f"Loading disruption file: {args.disruptions}")
    df = pd.read_csv(args.disruptions, sep="\t", low_memory=False)
    print(f"Rows: {len(df):,}")

    required = [
        "Chr", "Motif_Start", "Motif_End", "Name", "Gene",
        "Tumor_Sample_Barcode", "Variant_Chr", "Variant_Start",
        "Variant_End", "Reference_Allele", "Tumor_Seq_Allele2", "Entropy",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in disruption table: {missing}")

    gene_map = load_gene_map(args.gene_map)
    df = add_tf_and_gene_labels(df, gene_map)

    recurrent_loci, recurrent_loci_filtered = build_recurrent_loci(
        df,
        min_patients=args.min_patients,
        outdir=args.outdir,
    )

    gene_summary = build_gene_prioritization(
        df,
        recurrent_loci_filtered,
        args.outdir,
    )

    genes_of_interest = set(args.genes_of_interest) if args.genes_of_interest else DEFAULT_TOP_GENES_OF_INTEREST
    tfs_of_interest = set(args.tfs_of_interest) if args.tfs_of_interest else DEFAULT_TOP_TFS_OF_INTEREST

    integrated, gene_tf, tf_gene, intersection = build_integrated_locus_tables(
        df,
        recurrent_loci_filtered,
        args.outdir,
        top_n=args.top_n_items,
        genes_of_interest=genes_of_interest,
        tfs_of_interest=tfs_of_interest,
    )

    hotspots = build_hotspots(
        integrated,
        args.outdir,
        merge_gap=args.merge_gap,
        top_n_tfs=args.top_n_tfs,
        top_n_genes=args.top_n_genes,
    )

    motif_summary = build_hotspot_motif_summary(
        hotspots=hotspots,
        disruption_df=df,
        functional_bed=args.functional_bed,
        functional_table=args.functional_table,
        chunk_glob=args.chunk_glob,
        outdir=args.outdir,
        top_n_hotspots=args.top_hotspots,
    )

    non_tss = build_non_tss_examples(motif_summary, args.outdir)

    write_summary(
        outdir=args.outdir,
        df=df,
        recurrent_loci=recurrent_loci,
        recurrent_loci_filtered=recurrent_loci_filtered,
        gene_summary=gene_summary,
        integrated=integrated,
        hotspots=hotspots,
        motif_summary=motif_summary,
        non_tss=non_tss,
        min_patients=args.min_patients,
        merge_gap=args.merge_gap,
    )

    print("\nDone.")
    print(f"Outputs saved in: {args.outdir}")

    if not hotspots.empty:
        print("\nTop hotspots:")
        print(hotspots.head(10).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build recurrent hotspot, gene and TF prioritization tables."
    )

    parser.add_argument(
        "--disruptions",
        required=True,
        type=Path,
        help="Motif disrupting variants table from script 05."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory."
    )
    parser.add_argument(
        "--gene-map",
        type=Path,
        default=None,
        help="Optional gene ID to gene symbol mapping table."
    )
    parser.add_argument(
        "--min-patients",
        type=int,
        default=10,
        help="Minimum unique patients required for a locus to be considered recurrent."
    )
    parser.add_argument(
        "--merge-gap",
        type=int,
        default=200,
        help="Maximum gap in bp for merging recurrent loci into a hotspot."
    )
    parser.add_argument(
        "--top-n-items",
        type=int,
        default=5,
        help="Number of top genes/TFs to keep per recurrent locus."
    )
    parser.add_argument(
        "--top-n-tfs",
        type=int,
        default=5,
        help="Number of top TFs to report per hotspot."
    )
    parser.add_argument(
        "--top-n-genes",
        type=int,
        default=5,
        help="Number of top genes to report per hotspot."
    )
    parser.add_argument(
        "--top-hotspots",
        type=int,
        default=10,
        help="Number of top hotspots to use for optional motif summary."
    )
    parser.add_argument(
        "--functional-bed",
        type=Path,
        default=None,
        help="Optional functional motif BED for hotspot motif summaries."
    )
    parser.add_argument(
        "--functional-table",
        type=Path,
        default=None,
        help="Optional functional prediction table with mid and predictedClass columns."
    )
    parser.add_argument(
        "--chunk-glob",
        default=None,
        help="Optional glob for prepared ROSETTA input chunks to recover feature summaries."
    )
    parser.add_argument(
        "--genes-of-interest",
        nargs="*",
        default=None,
        help="Optional genes used for focused gene/TF intersection shortlist."
    )
    parser.add_argument(
        "--tfs-of-interest",
        nargs="*",
        default=None,
        help="Optional TFs used for focused gene/TF intersection shortlist."
    )

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
