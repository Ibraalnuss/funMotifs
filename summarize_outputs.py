#!/usr/bin/env python3

"""
summarize_outputs.py

Summarize the main outputs of the funMotifs CRC ROSETTA pipeline.

This script is useful for thesis traceability. It collects row counts and
important summary values from the main output tables and writes:

    results/tables/pipeline_summary.tsv
    results/tables/pipeline_summary.txt

Example:
    python scripts/summarize_outputs.py \
      --data-processed data/processed \
      --results-tables results/tables \
      --out-prefix results/tables/pipeline_summary
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def count_rows(path: Path, sep: str = "\t", has_header: bool = True) -> int | None:
    if not path.exists() or path.stat().st_size == 0:
        return None

    with path.open("r") as handle:
        n = sum(1 for _ in handle)

    if has_header and n > 0:
        return n - 1

    return n


def read_table(path: Path, sep: str = "\t") -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep=sep, low_memory=False)


def add_metric(records: list[dict], section: str, metric: str, value, source: Path | str = "") -> None:
    records.append(
        {
            "Section": section,
            "Metric": metric,
            "Value": "" if value is None else value,
            "Source": str(source),
        }
    )


def safe_nunique(df: pd.DataFrame, col: str):
    if df.empty or col not in df.columns:
        return None
    return df[col].nunique()


def safe_max(df: pd.DataFrame, col: str):
    if df.empty or col not in df.columns:
        return None
    value = pd.to_numeric(df[col], errors="coerce").max()
    if pd.isna(value):
        return None
    return value


def summarize(args: argparse.Namespace) -> pd.DataFrame:
    records: list[dict] = []

    data_processed = args.data_processed
    tables = args.results_tables

    # ------------------------------------------------------------------
    # ROSETTA prediction outputs
    # ------------------------------------------------------------------
    functional_table = data_processed / "colon_rosetta_functional_only.tsv"
    functional_bed = data_processed / "colon_rosetta_functional_motifs.sorted.bed"
    full_predictions = data_processed / "colon_rosetta_predictions.tsv"
    merge_summary = tables / "rosetta_merge_summary.tsv"

    add_metric(records, "ROSETTA prediction", "functional motif table rows", count_rows(functional_table), functional_table)
    add_metric(records, "ROSETTA prediction", "functional motif BED rows", count_rows(functional_bed, has_header=False), functional_bed)
    add_metric(records, "ROSETTA prediction", "full prediction rows", count_rows(full_predictions), full_predictions)

    merge = read_table(merge_summary)
    if not merge.empty:
        for col in ["functional_files", "full_prediction_files", "bed_files", "functional_rows", "full_prediction_rows", "sorted_bed_rows"]:
            if col in merge.columns:
                add_metric(records, "ROSETTA prediction", col, merge[col].iloc[0], merge_summary)

    # ------------------------------------------------------------------
    # Variant overlap and disruption
    # ------------------------------------------------------------------
    variant_dir = tables / "variant_overlap_disruption"
    overlap_file = variant_dir / "functional_motif_variant_overlaps.tsv"
    disruption_file = variant_dir / "motif_disrupting_variants.tsv"
    cleaned_dir = variant_dir / "cleaned"

    unique_patient_file = cleaned_dir / "unique_motif_region_per_patient.tsv"
    unique_loci_file = cleaned_dir / "collapsed_disrupted_loci.tsv"
    funnel_file = cleaned_dir / "overlap_disruption_funnel.tsv"

    add_metric(records, "Variant overlap", "motif variant overlap rows", count_rows(overlap_file), overlap_file)
    add_metric(records, "Variant overlap", "retained motif disrupting rows", count_rows(disruption_file), disruption_file)
    add_metric(records, "Variant overlap", "unique motif region per patient rows", count_rows(unique_patient_file), unique_patient_file)
    add_metric(records, "Variant overlap", "collapsed disrupted loci rows", count_rows(unique_loci_file), unique_loci_file)

    disruptions = read_table(disruption_file)
    add_metric(records, "Variant overlap", "unique patients with retained disruptions", safe_nunique(disruptions, "Tumor_Sample_Barcode"), disruption_file)
    add_metric(records, "Variant overlap", "unique disrupted TF motif names", safe_nunique(disruptions, "Name"), disruption_file)
    add_metric(records, "Variant overlap", "maximum disruption score or entropy", safe_max(disruptions, "Entropy"), disruption_file)

    funnel = read_table(funnel_file)
    if not funnel.empty and {"Step", "Count"}.issubset(funnel.columns):
        for _, row in funnel.iterrows():
            add_metric(records, "Variant overlap funnel", str(row["Step"]), row["Count"], funnel_file)

    # ------------------------------------------------------------------
    # Recurrent hotspots and genes
    # ------------------------------------------------------------------
    recurrent_dir = tables / "recurrent_hotspots_and_genes"
    recurrent_loci = recurrent_dir / "recurrent_loci_filtered.tsv"
    hotspots_file = recurrent_dir / "recurrent_hotspot_regions.tsv"
    gene_table_file = recurrent_dir / "gene_prioritization_from_recurrent_loci.tsv"
    integrated_file = recurrent_dir / "integrated_recurrent_locus_table.tsv"
    non_tss_file = recurrent_dir / "non_tss_hotspot_examples.tsv"

    add_metric(records, "Recurrent analysis", "recurrent loci rows", count_rows(recurrent_loci), recurrent_loci)
    add_metric(records, "Recurrent analysis", "recurrent hotspot rows", count_rows(hotspots_file), hotspots_file)
    add_metric(records, "Recurrent analysis", "prioritized gene rows", count_rows(gene_table_file), gene_table_file)
    add_metric(records, "Recurrent analysis", "integrated locus rows", count_rows(integrated_file), integrated_file)
    add_metric(records, "Recurrent analysis", "non TSS hotspot example rows", count_rows(non_tss_file), non_tss_file)

    hotspots = read_table(hotspots_file)
    add_metric(records, "Recurrent analysis", "max hotspot unique patients", safe_max(hotspots, "Max_Unique_Patients"), hotspots_file)
    add_metric(records, "Recurrent analysis", "max loci per hotspot", safe_max(hotspots, "Num_Loci_In_Hotspot"), hotspots_file)

    genes = read_table(gene_table_file)
    add_metric(records, "Recurrent analysis", "unique prioritized genes", safe_nunique(genes, "Gene_Label") or safe_nunique(genes, "Gene"), gene_table_file)
    add_metric(records, "Recurrent analysis", "max gene unique patients", safe_max(genes, "Unique_Patients"), gene_table_file)

    # ------------------------------------------------------------------
    # External interpretation
    # ------------------------------------------------------------------
    external_dir = tables / "external_interpretation"
    gwas_besthit = external_dir / "crc_gwas_hotspot_besthit.tsv"
    gwas_shortlist = external_dir / "crc_gwas_hotspot_besthit_shortlist.tsv"
    exact_gwas = external_dir / "exact_variant_gwas_matches.tsv"
    exact_gwas_collapsed = external_dir / "exact_variant_gwas_matches_collapsed.tsv"
    ncg_overlap = external_dir / "top500_genes_overlapping_NCG.tsv"
    top_pathway = external_dir / "top_genes_pathway_table.tsv"
    ncg_pathway = external_dir / "ncg_genes_pathway_table.tsv"

    add_metric(records, "External interpretation", "GWAS besthit rows", count_rows(gwas_besthit), gwas_besthit)
    add_metric(records, "External interpretation", "GWAS shortlist rows", count_rows(gwas_shortlist), gwas_shortlist)
    add_metric(records, "External interpretation", "exact variant GWAS match rows", count_rows(exact_gwas), exact_gwas)
    add_metric(records, "External interpretation", "collapsed exact variant GWAS match rows", count_rows(exact_gwas_collapsed), exact_gwas_collapsed)
    add_metric(records, "External interpretation", "NCG overlap rows", count_rows(ncg_overlap), ncg_overlap)
    add_metric(records, "External interpretation", "top gene pathway rows", count_rows(top_pathway), top_pathway)
    add_metric(records, "External interpretation", "NCG gene pathway rows", count_rows(ncg_pathway), ncg_pathway)

    gwas = read_table(gwas_besthit)
    add_metric(records, "External interpretation", "GWAS exact or <=1kb hotspots", count_gwas_tier(gwas, ["Tier 1: exact", "Tier 1: <=1kb"]), gwas_besthit)
    add_metric(records, "External interpretation", "GWAS gene concordant hotspots", count_gene_concordance(gwas), gwas_besthit)

    # ------------------------------------------------------------------
    # Patient subgroup analysis
    # ------------------------------------------------------------------
    subgroup_dir = tables / "patient_hotspot_subgroups"
    matrix_file = subgroup_dir / "patient_by_hotspot_binary_matrix.tsv"
    cluster_file = subgroup_dir / "patient_jaccard_clusters.tsv"
    cluster_sizes_file = subgroup_dir / "patient_cluster_sizes.tsv"
    profile_file = subgroup_dir / "cluster_hotspot_profiles.tsv"
    gene_summary_file = subgroup_dir / "cluster_gene_tf_summary.tsv"

    add_metric(records, "Patient subgroups", "patient by hotspot matrix rows", count_rows(matrix_file), matrix_file)
    add_metric(records, "Patient subgroups", "patient cluster assignment rows", count_rows(cluster_file), cluster_file)
    add_metric(records, "Patient subgroups", "cluster size rows", count_rows(cluster_sizes_file), cluster_sizes_file)
    add_metric(records, "Patient subgroups", "cluster hotspot profile rows", count_rows(profile_file), profile_file)
    add_metric(records, "Patient subgroups", "cluster gene TF summary rows", count_rows(gene_summary_file), gene_summary_file)

    clusters = read_table(cluster_file)
    add_metric(records, "Patient subgroups", "number of patient subgroups", safe_nunique(clusters, "Cluster_Label"), cluster_file)

    if not clusters.empty and "Jaccard_Threshold" in clusters.columns:
        add_metric(records, "Patient subgroups", "chosen Jaccard threshold", clusters["Jaccard_Threshold"].iloc[0], cluster_file)

    cluster_sizes = read_table(cluster_sizes_file)
    add_metric(records, "Patient subgroups", "largest subgroup size", safe_max(cluster_sizes, "Num_Patients"), cluster_sizes_file)

    return pd.DataFrame(records)


def count_gwas_tier(df: pd.DataFrame, tiers: list[str]):
    if df.empty or "Evidence_Tier" not in df.columns:
        return None
    return int(df["Evidence_Tier"].isin(tiers).sum())


def count_gene_concordance(df: pd.DataFrame):
    if df.empty or "Gene_Concordance" not in df.columns:
        return None
    return int(df["Gene_Concordance"].astype(str).str.lower().eq("yes").sum())


def write_text_summary(summary: pd.DataFrame, output: Path) -> None:
    lines = ["Pipeline output summary", ""]

    for section, sub in summary.groupby("Section", sort=False):
        lines.append(section)
        lines.append("-" * len(section))

        for _, row in sub.iterrows():
            value = row["Value"]
            metric = row["Metric"]
            source = row["Source"]
            lines.append(f"{metric}: {value}  [{source}]")

        lines.append("")

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize main outputs from the pipeline."
    )

    parser.add_argument(
        "--data-processed",
        type=Path,
        default=Path("data/processed"),
        help="Processed data directory."
    )
    parser.add_argument(
        "--results-tables",
        type=Path,
        default=Path("results/tables"),
        help="Results table directory."
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("results/tables/pipeline_summary"),
        help="Output prefix without extension."
    )

    args = parser.parse_args()

    summary = summarize(args)

    out_tsv = args.out_prefix.with_suffix(".tsv")
    out_txt = args.out_prefix.with_suffix(".txt")

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_tsv, sep="\t", index=False)
    write_text_summary(summary, out_txt)

    print(f"Wrote: {out_tsv}")
    print(f"Wrote: {out_txt}")

    print("\nSummary preview:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
