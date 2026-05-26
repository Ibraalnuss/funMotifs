#!/usr/bin/env python3

"""
05_variant_overlap_and_disruption.py

Overlap predicted functional colon motifs with CRC somatic variants and
identify candidate motif-disrupting events.

Steps:
  1. Sort the functional motif BED and (optionally) the variant BED by coordinate.
  2. Run bedtools intersect to find all motif-variant overlaps.
  3. For SNPs: retain only those where the reference allele frequency minus the
     alternate allele frequency (from the motif PFM) exceeds --entropy-threshold.
     Non-SNP variants (indels, DNPs, etc.) are kept automatically.
  4. Deduplicate and summarize: per-patient unique motif regions, per-locus
     recurrence counts, top TFs by patient count, multi-motif clusters.

Requires:
    bedtools in PATH
    A motif PFM table with columns: name, position, allele, freq

Example:
    python scripts/05_variant_overlap_and_disruption.py \\
      --functional-bed data/processed/colon_rosetta_functional_motifs.sorted.bed \\
      --variant-bed data/external/CRC-colon.section6.sorted.bed \\
      --pfm data/external/motifs_pfm.tsv \\
      --outdir results/tables/variant_overlap_disruption
"""

from __future__ import annotations

import argparse
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


OVERLAP_COLUMNS = [
    "Chr",
    "Motif_Start",
    "Motif_End",
    "Name",
    "Strand",
    "mid",
    "Variant_Chr",
    "Variant_Start",
    "Variant_End",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
    "Transcript_ID",
    "Gene",
    "Overlapping_Base_Pairs",
]


def check_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    if path.stat().st_size == 0:
        raise RuntimeError(f"{label} exists but is empty: {path}")


def sort_bed(input_bed: Path, output_bed: Path) -> None:
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    with output_bed.open("w") as out:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", str(input_bed)],
            check=True,
            stdout=out,
        )

    check_file(output_bed, "sorted BED")
    print(f"Saved sorted BED: {output_bed}")


def run_bedtools_intersect(motif_bed: Path, variant_bed: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w") as out:
        subprocess.run(
            [
                "bedtools",
                "intersect",
                "-sorted",
                "-a",
                str(motif_bed),
                "-b",
                str(variant_bed),
                "-wo",
            ],
            check=True,
            stdout=out,
        )

    check_file(output, "motif-variant overlap file")
    print(f"Saved overlap file: {output}")


def validate_overlap_file(path: Path) -> None:
    n_rows = 0
    bad_rows = 0

    with path.open("r") as handle:
        for line in handle:
            n_rows += 1
            if len(line.rstrip("\n").split("\t")) != len(OVERLAP_COLUMNS):
                bad_rows += 1

    print(f"Overlap rows: {n_rows:,}")
    print(f"Rows with wrong column count: {bad_rows:,}")

    if bad_rows:
        raise RuntimeError(
            f"Overlap file has {bad_rows:,} rows with wrong column count. "
            f"Expected {len(OVERLAP_COLUMNS)} columns."
        )


def compute_position_in_motif(row: pd.Series) -> int:
    """
    Compute the variant position relative to the motif.

    Variant_Start is BED start, meaning 0-based.
    """
    if row["Strand"] == "+":
        return int(row["Variant_Start"]) - int(row["Motif_Start"])
    return int(row["Motif_End"]) - int(row["Variant_Start"]) - 1


def load_pfm_lookup(pfm_path: Path) -> dict[tuple[str, int, str], float]:
    pfm = pd.read_csv(pfm_path, sep="\t")

    required = {"name", "position", "allele", "freq"}
    missing = required - set(pfm.columns)
    if missing:
        raise ValueError(f"PFM file is missing columns: {sorted(missing)}")

    pfm["name"] = pfm["name"].astype(str)
    pfm["position"] = pfm["position"].astype(int)
    pfm["allele"] = pfm["allele"].astype(str).str.upper()

    return {
        (row["name"], row["position"], row["allele"]): float(row["freq"])
        for _, row in pfm.iterrows()
    }


def call_disrupting_variants(
    overlap_path: Path,
    pfm_path: Path,
    output_path: Path,
    entropy_threshold: float,
) -> pd.DataFrame:
    print("Loading overlap file...")
    overlap = pd.read_csv(
        overlap_path,
        sep="\t",
        header=None,
        names=OVERLAP_COLUMNS,
        low_memory=False,
    )

    print("Loading PFM lookup...")
    pfm_lookup = load_pfm_lookup(pfm_path)

    rows = []
    n_non_snp = 0
    n_snp = 0
    n_missing_lookup = 0
    n_kept_snp = 0

    for _, row in overlap.iterrows():
        row = row.copy()
        row["Reference_Allele"] = str(row["Reference_Allele"]).upper()
        row["Tumor_Seq_Allele2"] = str(row["Tumor_Seq_Allele2"]).upper()

        if row["Variant_Type"] != "SNP":
            row["Entropy"] = 1.0
            rows.append(row)
            n_non_snp += 1
            continue

        n_snp += 1
        motif_position = compute_position_in_motif(row)

        ref_key = (row["Name"], motif_position, row["Reference_Allele"])
        alt_key = (row["Name"], motif_position, row["Tumor_Seq_Allele2"])

        ref_freq = pfm_lookup.get(ref_key)
        alt_freq = pfm_lookup.get(alt_key)

        if ref_freq is None or alt_freq is None:
            n_missing_lookup += 1
            continue

        row["Entropy"] = ref_freq - alt_freq

        if row["Entropy"] > entropy_threshold:
            rows.append(row)
            n_kept_snp += 1

    out = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False)

    print(f"Saved disrupting variants: {output_path}")
    print(f"Input overlap rows: {len(overlap):,}")
    print(f"Non-SNP rows kept automatically: {n_non_snp:,}")
    print(f"SNP rows evaluated: {n_snp:,}")
    print(f"SNP rows missing PFM lookup: {n_missing_lookup:,}")
    print(f"SNP rows kept by entropy threshold: {n_kept_snp:,}")
    print(f"Final disrupting rows: {len(out):,}")

    return out


def save_df(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def remove_duplicate_rows(df: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    subset_cols = [col for col in df.columns if col != "mid"]
    out = df.drop_duplicates(subset=subset_cols).copy()
    save_df(out, outdir / "motif_disrupting_variants_unique.tsv")
    return out


def make_unique_motif_region_per_patient(df: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    subset_cols = ["Chr", "Motif_Start", "Motif_End", "Name", "Tumor_Sample_Barcode"]
    out = df.drop_duplicates(subset=subset_cols).copy()
    save_df(out, outdir / "unique_motif_region_per_patient.tsv")
    return out


def make_unique_motif_regions(df: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    subset_cols = ["Chr", "Motif_Start", "Motif_End", "Name", "Transcript_ID", "Gene"]
    out = df.drop_duplicates(subset=subset_cols).copy()
    save_df(out, outdir / "unique_motif_regions.tsv")
    return out


def collapse_disrupted_loci(df_patient: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    group_cols = ["Chr", "Motif_Start", "Motif_End"]

    count_df = (
        df_patient.groupby(group_cols)
        .size()
        .reset_index(name="Num_MotifRegionPatient_Records")
    )

    agg_df = (
        df_patient.groupby(group_cols)
        .agg(lambda x: ",".join(map(str, sorted(set(map(str, x))))))
        .reset_index()
    )

    out = agg_df.merge(count_df, on=group_cols)
    save_df(out, outdir / "collapsed_disrupted_loci.tsv")
    return out


def summarize_top_tfs(df_patient: pd.DataFrame, outdir: Path, top_n: int) -> pd.DataFrame:
    tf_patients = defaultdict(set)

    for _, row in df_patient.iterrows():
        motif_name = str(row["Name"])
        tf = motif_name.split("_")[0]
        patient = str(row["Tumor_Sample_Barcode"])
        tf_patients[tf].add(patient)

    out = pd.DataFrame(
        [{"TF": tf, "Unique_Patients": len(patients)} for tf, patients in tf_patients.items()]
    )

    if not out.empty:
        out = out.sort_values("Unique_Patients", ascending=False).head(top_n)

    save_df(out, outdir / "top_tfs_by_unique_patients.tsv")
    return out


def summarize_top_motif_regions(df_patient: pd.DataFrame, outdir: Path, top_n: int) -> pd.DataFrame:
    region_patients = defaultdict(set)

    for _, row in df_patient.iterrows():
        region = (str(row["Chr"]), int(row["Motif_Start"]), int(row["Motif_End"]))
        patient = str(row["Tumor_Sample_Barcode"])
        region_patients[region].add(patient)

    out = pd.DataFrame(
        [
            {
                "Chr": region[0],
                "Motif_Start": region[1],
                "Motif_End": region[2],
                "Unique_Patients": len(patients),
            }
            for region, patients in region_patients.items()
        ]
    )

    if not out.empty:
        out = out.sort_values("Unique_Patients", ascending=False).head(top_n)

    save_df(out, outdir / "top_motif_regions_by_unique_patients.tsv")
    return out


def cluster_motif_regions(df_regions: pd.DataFrame, outdir: Path, distance: int) -> pd.DataFrame:
    motifs = []

    for _, row in df_regions.iterrows():
        chrom = str(row["Chr"])
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom

        motifs.append(
            {
                "chr": chrom,
                "start": int(row["Motif_Start"]),
                "end": int(row["Motif_End"]),
                "id": str(row["Name"]),
                "transcript_id": str(row["Transcript_ID"]),
                "gene": str(row["Gene"]),
            }
        )

    motifs.sort(key=lambda x: (x["chr"], x["start"], x["end"], x["id"]))

    clusters = []
    current = []

    for motif in motifs:
        if not current:
            current.append(motif)
            continue

        last = current[-1]
        if motif["chr"] == last["chr"] and (motif["start"] - last["end"]) <= distance:
            current.append(motif)
        else:
            clusters.append(current)
            current = [motif]

    if current:
        clusters.append(current)

    rows = []
    cluster_id = 0

    for cluster in clusters:
        if len(cluster) <= 1:
            continue

        cluster_id += 1
        cluster_chr = cluster[0]["chr"]
        cluster_start = min(m["start"] for m in cluster)
        cluster_end = max(m["end"] for m in cluster)

        for motif in cluster:
            rows.append(
                {
                    "Cluster_ID": f"Cluster_{cluster_id}",
                    "Cluster_Chr": cluster_chr,
                    "Cluster_Start": cluster_start,
                    "Cluster_End": cluster_end,
                    "Num_Motifs_in_Cluster": len(cluster),
                    "Motif_ID": motif["id"],
                    "Motif_Chr": motif["chr"],
                    "Motif_Start": motif["start"],
                    "Motif_End": motif["end"],
                    "Transcript_ID": motif["transcript_id"],
                    "Gene": motif["gene"],
                }
            )

    out = pd.DataFrame(rows)
    save_df(out, outdir / "multi_motif_clusters.tsv")
    return out


def write_analysis_summary(
    outdir: Path,
    overlap_rows: int,
    disruption_rows: int,
    unique_rows: pd.DataFrame,
    patient_rows: pd.DataFrame,
    region_rows: pd.DataFrame,
    collapsed_rows: pd.DataFrame,
    top_tfs: pd.DataFrame,
    top_regions: pd.DataFrame,
    clusters: pd.DataFrame,
) -> None:
    summary = pd.DataFrame(
        [
            ["functional_motif_variant_overlaps", overlap_rows],
            ["retained_disrupting_rows", disruption_rows],
            ["unique_rows_excluding_mid", len(unique_rows)],
            ["unique_motif_region_per_patient", len(patient_rows)],
            ["unique_motif_regions", len(region_rows)],
            ["unique_disrupted_loci", len(collapsed_rows)],
            [
                "top_tf_max_unique_patients",
                int(top_tfs["Unique_Patients"].max()) if not top_tfs.empty else 0,
            ],
            [
                "top_motif_region_max_unique_patients",
                int(top_regions["Unique_Patients"].max()) if not top_regions.empty else 0,
            ],
            [
                "clusters_with_more_than_one_motif",
                int(clusters["Cluster_ID"].nunique()) if not clusters.empty else 0,
            ],
        ],
        columns=["Metric", "Value"],
    )

    save_df(summary, outdir / "analysis_summary.tsv")


def write_funnel(outdir: Path, summary_path: Path) -> None:
    summary = pd.read_csv(summary_path, sep="\t")
    summary_map = dict(zip(summary["Metric"], summary["Value"]))

    rows = [
        ["functional_motif_variant_overlaps", int(summary_map.get("functional_motif_variant_overlaps", 0))],
        ["retained_disrupting_rows", int(summary_map.get("retained_disrupting_rows", 0))],
        ["unique_rows_excluding_mid", int(summary_map.get("unique_rows_excluding_mid", 0))],
        ["unique_motif_region_per_patient", int(summary_map.get("unique_motif_region_per_patient", 0))],
        ["unique_disrupted_loci", int(summary_map.get("unique_disrupted_loci", 0))],
    ]

    funnel = pd.DataFrame(rows, columns=["Step", "Count"])
    funnel["Percent_of_previous"] = funnel["Count"] / funnel["Count"].shift(1) * 100
    funnel.loc[0, "Percent_of_previous"] = 100.0
    funnel["Percent_of_start"] = funnel["Count"] / funnel.loc[0, "Count"] * 100

    save_df(funnel, outdir / "overlap_disruption_funnel.tsv")


def run_analysis(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)

    sorted_motifs = args.outdir / "functional_motifs.sorted.bed"
    overlap_path = args.outdir / "functional_motif_variant_overlaps.tsv"
    disruption_path = args.outdir / "motif_disrupting_variants.tsv"

    check_file(args.functional_bed, "functional motif BED")
    check_file(args.variant_bed, "variant BED")
    check_file(args.pfm, "PFM table")

    sort_bed(args.functional_bed, sorted_motifs)

    if args.variant_bed_is_sorted:
        sorted_variants = args.variant_bed
    else:
        sorted_variants = args.outdir / "variants.sorted.bed"
        sort_bed(args.variant_bed, sorted_variants)

    run_bedtools_intersect(sorted_motifs, sorted_variants, overlap_path)
    validate_overlap_file(overlap_path)

    disruptions = call_disrupting_variants(
        overlap_path=overlap_path,
        pfm_path=args.pfm,
        output_path=disruption_path,
        entropy_threshold=args.entropy_threshold,
    )

    print("\nCleaning and summarizing disruption events...")
    clean_dir = args.outdir / "cleaned"
    clean_dir.mkdir(parents=True, exist_ok=True)

    unique = remove_duplicate_rows(disruptions, clean_dir)
    patient = make_unique_motif_region_per_patient(unique, clean_dir)
    regions = make_unique_motif_regions(unique, clean_dir)
    collapsed = collapse_disrupted_loci(patient, clean_dir)
    top_tfs = summarize_top_tfs(patient, clean_dir, args.top_tfs)
    top_regions = summarize_top_motif_regions(patient, clean_dir, args.top_regions)
    clusters = cluster_motif_regions(regions, clean_dir, args.cluster_distance)

    summary_path = clean_dir / "analysis_summary.tsv"

    write_analysis_summary(
        outdir=clean_dir,
        overlap_rows=sum(1 for _ in overlap_path.open()),
        disruption_rows=len(disruptions),
        unique_rows=unique,
        patient_rows=patient,
        region_rows=regions,
        collapsed_rows=collapsed,
        top_tfs=top_tfs,
        top_regions=top_regions,
        clusters=clusters,
    )

    write_funnel(clean_dir, summary_path)

    print("\nDone.")
    print(f"Overlap file: {overlap_path}")
    print(f"Disruption file: {disruption_path}")
    print(f"Cleaned outputs: {clean_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Overlap functional motifs with CRC variants and call motif-disrupting events."
    )

    parser.add_argument(
        "--functional-bed",
        required=True,
        type=Path,
        help="Functional motif BED file from ROSETTA prediction."
    )
    parser.add_argument(
        "--variant-bed",
        required=True,
        type=Path,
        help="CRC variant BED file. Must contain the 10 expected variant columns."
    )
    parser.add_argument(
        "--pfm",
        required=True,
        type=Path,
        help="Motif PFM table with columns: name, position, allele, freq."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory."
    )
    parser.add_argument(
        "--variant-bed-is-sorted",
        action="store_true",
        help="Use this if the variant BED is already sorted."
    )
    parser.add_argument(
        "--entropy-threshold",
        type=float,
        default=0.3,
        help="Minimum ref minus alt PFM frequency difference for SNP retention."
    )
    parser.add_argument(
        "--top-tfs",
        type=int,
        default=20,
        help="Number of top TFs to report."
    )
    parser.add_argument(
        "--top-regions",
        type=int,
        default=10,
        help="Number of top motif regions to report."
    )
    parser.add_argument(
        "--cluster-distance",
        type=int,
        default=200,
        help="Maximum distance between motif regions to place them in the same cluster."
    )

    args = parser.parse_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
