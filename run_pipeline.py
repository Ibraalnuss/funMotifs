#!/usr/bin/env python3

"""
run_pipeline.py

Convenience wrapper that runs the numbered pipeline scripts in order.

This script contains no analysis logic itself -- it just calls the right
script with the right arguments and checks that each step succeeds before
moving to the next.

Run a single step:
    python scripts/run_pipeline.py --step rosetta_prepare

Run a range of steps:
    python scripts/run_pipeline.py --from rosetta_prepare --to external

Full run (requires all external data files to be in place):
    python scripts/run_pipeline.py --from rosetta_prepare --to subgroups \\
      --colon-annotations data/external/colon_annotations.tsv \\
      --rules data/external/sig_rules_final_with_pretty_rules.rds \\
      --deployment data/external/rosetta_deployment_info.rds \\
      --variant-bed data/external/CRC-colon.section6.sorted.bed \\
      --pfm data/external/motifs_pfm.tsv \\
      --gwas data/external/gwas-association-CRC.tsv \\
      --ncg data/external/NCG_cancerdrivers_annotation_supporting_evidence.tsv

Use --dry-run to print all commands without running them.

Pipeline step order:
    decision_table -> rosetta_prepare -> rosetta_predict -> rosetta_merge ->
    variant_overlap -> recurrent -> external -> figures -> subgroups
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent

STEP_ORDER = [
    "decision_table",
    "rosetta_prepare",
    "rosetta_predict",
    "rosetta_merge",
    "variant_overlap",
    "recurrent",
    "external",
    "figures",
    "subgroups",
]


def run_command(cmd: list[str], dry_run: bool = False) -> None:
    print("\n" + "=" * 80)
    print("Running:")
    print(" ".join(str(x) for x in cmd))
    print("=" * 80)

    if dry_run:
        return

    subprocess.run(cmd, check=True)


def require_path(path: Path | None, label: str) -> Path:
    if path is None:
        raise ValueError(f"Missing required argument: {label}")
    return path


def selected_steps(args: argparse.Namespace) -> list[str]:
    if args.step:
        return [args.step]

    start = STEP_ORDER.index(args.from_step)
    end = STEP_ORDER.index(args.to_step)

    if start > end:
        raise ValueError("--from must come before --to in the pipeline order.")

    return STEP_ORDER[start : end + 1]


def build_commands(args: argparse.Namespace) -> dict[str, list[str]]:
    python = sys.executable

    paths = {
        "colon_input": args.data_interim / "colon_rosetta_input.tsv",
        "chunk_dir": args.data_interim / "colon_rosetta_chunks",
        "chunk_list": args.data_interim / "colon_rosetta_chunk_list.txt",
        "rosetta_predictions": args.results_tables / "rosetta_chunk_predictions",
        "functional_tsv": args.data_processed / "colon_rosetta_functional_only.tsv",
        "full_predictions_tsv": args.data_processed / "colon_rosetta_predictions.tsv",
        "functional_bed": args.data_processed / "colon_rosetta_functional_motifs.bed",
        "functional_sorted_bed": args.data_processed / "colon_rosetta_functional_motifs.sorted.bed",
        "merge_summary": args.results_tables / "rosetta_merge_summary.tsv",
        "variant_outdir": args.results_tables / "variant_overlap_disruption",
        "recurrent_outdir": args.results_tables / "recurrent_hotspots_and_genes",
        "external_outdir": args.results_tables / "external_interpretation",
        "figures_outdir": args.results_figures,
        "subgroup_outdir": args.results_tables / "patient_hotspot_subgroups",
        "subgroup_figdir": args.results_figures / "patient_hotspot_subgroups",
    }

    commands: dict[str, list[str]] = {}

    commands["decision_table"] = [
        python,
        str(SCRIPT_DIR / "01_prepare_decision_table.py"),
        "evaluate",
        "--table",
        str(require_path(args.decision_table, "--decision-table")),
        "--outdir",
        str(args.results_tables / "decision_table_evaluation"),
        "--name",
        args.decision_table_name,
    ]

    commands["rosetta_prepare"] = [
        python,
        str(SCRIPT_DIR / "02_prepare_rosetta_colon_input.py"),
        "--input",
        str(require_path(args.colon_annotations, "--colon-annotations")),
        "--output",
        str(paths["colon_input"]),
        "--chunk-dir",
        str(paths["chunk_dir"]),
        "--chunk-size",
        str(args.chunk_size),
        "--chunk-list",
        str(paths["chunk_list"]),
    ]

    commands["rosetta_predict"] = [
        "Rscript",
        str(SCRIPT_DIR / "03_predict_rosetta_colon.R"),
        "--chunk-list",
        str(paths["chunk_list"]),
        "--rules",
        str(require_path(args.rules, "--rules")),
        "--deployment",
        str(require_path(args.deployment, "--deployment")),
        "--out-dir",
        str(paths["rosetta_predictions"]),
    ]

    commands["rosetta_merge"] = [
        python,
        str(SCRIPT_DIR / "04_merge_rosetta_predictions.py"),
        "--input-dir",
        str(paths["rosetta_predictions"]),
        "--out-functional",
        str(paths["functional_tsv"]),
        "--out-full",
        str(paths["full_predictions_tsv"]),
        "--out-bed",
        str(paths["functional_bed"]),
        "--out-sorted-bed",
        str(paths["functional_sorted_bed"]),
        "--out-summary",
        str(paths["merge_summary"]),
    ]

    commands["variant_overlap"] = [
        python,
        str(SCRIPT_DIR / "05_variant_overlap_and_disruption.py"),
        "--functional-bed",
        str(paths["functional_sorted_bed"]),
        "--variant-bed",
        str(require_path(args.variant_bed, "--variant-bed")),
        "--pfm",
        str(require_path(args.pfm, "--pfm")),
        "--outdir",
        str(paths["variant_outdir"]),
        "--entropy-threshold",
        str(args.entropy_threshold),
    ]

    if args.variant_bed_is_sorted:
        commands["variant_overlap"].append("--variant-bed-is-sorted")

    commands["recurrent"] = [
        python,
        str(SCRIPT_DIR / "06_recurrent_hotspots_and_genes.py"),
        "--disruptions",
        str(paths["variant_outdir"] / "motif_disrupting_variants.tsv"),
        "--outdir",
        str(paths["recurrent_outdir"]),
        "--min-patients",
        str(args.min_recurrent_patients),
        "--merge-gap",
        str(args.hotspot_merge_gap),
        "--functional-bed",
        str(paths["functional_sorted_bed"]),
        "--functional-table",
        str(paths["functional_tsv"]),
        "--chunk-glob",
        str(paths["chunk_dir"] / "colon_rosetta_input_chunk_*.tsv"),
    ]

    if args.gene_map:
        commands["recurrent"].extend(["--gene-map", str(args.gene_map)])

    commands["external"] = [
        python,
        str(SCRIPT_DIR / "07_external_interpretation.py"),
        "--hotspots",
        str(paths["recurrent_outdir"] / "recurrent_hotspot_regions.tsv"),
        "--genes",
        str(paths["recurrent_outdir"] / "gene_prioritization_from_recurrent_loci.tsv"),
        "--disruptions",
        str(paths["variant_outdir"] / "motif_disrupting_variants.tsv"),
        "--gwas",
        str(require_path(args.gwas, "--gwas")),
        "--outdir",
        str(paths["external_outdir"]),
    ]

    if args.ncg:
        commands["external"].extend(["--ncg", str(args.ncg)])

    if args.ld_overlap:
        commands["external"].extend(["--ld-overlap", str(args.ld_overlap)])

    if args.enrichr_top_genes:
        commands["external"].extend(["--enrichr-top-genes", str(args.enrichr_top_genes)])

    if args.enrichr_ncg_genes:
        commands["external"].extend(["--enrichr-ncg-genes", str(args.enrichr_ncg_genes)])

    commands["figures"] = [
        python,
        str(SCRIPT_DIR / "08_make_thesis_figures.py"),
        "--group",
        args.figure_group,
        "--recurrent-dir",
        str(paths["recurrent_outdir"]),
        "--external-dir",
        str(paths["external_outdir"]),
        "--variant-dir",
        str(paths["variant_outdir"]),
        "--functional-table",
        str(paths["functional_tsv"]),
        "--outdir",
        str(paths["figures_outdir"]),
    ]

    commands["subgroups"] = [
        python,
        str(SCRIPT_DIR / "09_patient_hotspot_subgroups.py"),
        "--hotspot-table",
        str(paths["recurrent_outdir"] / "recurrent_hotspot_regions.tsv"),
        "--event-table",
        str(paths["variant_outdir"] / "cleaned" / "unique_motif_region_per_patient.tsv"),
        "--outdir",
        str(paths["subgroup_outdir"]),
        "--figdir",
        str(paths["subgroup_figdir"]),
        "--min-patients",
        str(args.subgroup_min_patients),
        "--thresholds",
        args.subgroup_thresholds,
        "--target-main-groups",
        str(args.subgroup_target_main_groups),
        "--min-main-group-size",
        str(args.subgroup_min_main_group_size),
        "--min-main-group-fraction",
        str(args.subgroup_min_main_group_fraction),
        "--top-main-clusters",
        str(args.subgroup_top_main_clusters),
        "--top-hotspots-per-group",
        str(args.subgroup_top_hotspots_per_group),
        "--signature-threshold",
        str(args.subgroup_signature_threshold),
    ]

    if args.blacklist_bed:
        commands["subgroups"].extend(["--blacklist-bed", str(args.blacklist_bed)])

    return commands


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the final funMotifs CRC ROSETTA pipeline."
    )

    parser.add_argument(
        "--step",
        choices=STEP_ORDER,
        default=None,
        help="Run a single pipeline step."
    )
    parser.add_argument(
        "--from",
        dest="from_step",
        choices=STEP_ORDER,
        default="rosetta_prepare",
        help="First step to run when --step is not used."
    )
    parser.add_argument(
        "--to",
        dest="to_step",
        choices=STEP_ORDER,
        default="subgroups",
        help="Last step to run when --step is not used."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running them."
    )

    parser.add_argument("--data-external", type=Path, default=Path("data/external"))
    parser.add_argument("--data-interim", type=Path, default=Path("data/interim"))
    parser.add_argument("--data-processed", type=Path, default=Path("data/processed"))
    parser.add_argument("--results-tables", type=Path, default=Path("results/tables"))
    parser.add_argument("--results-figures", type=Path, default=Path("results/figures"))

    parser.add_argument("--decision-table", type=Path, default=None)
    parser.add_argument("--decision-table-name", default="decision_table")

    parser.add_argument("--colon-annotations", type=Path, default=None)
    parser.add_argument("--rules", type=Path, default=None)
    parser.add_argument("--deployment", type=Path, default=None)
    parser.add_argument("--variant-bed", type=Path, default=None)
    parser.add_argument("--variant-bed-is-sorted", action="store_true")
    parser.add_argument("--pfm", type=Path, default=None)
    parser.add_argument("--gene-map", type=Path, default=None)
    parser.add_argument("--gwas", type=Path, default=None)
    parser.add_argument("--ncg", type=Path, default=None)
    parser.add_argument("--ld-overlap", type=Path, default=None)
    parser.add_argument("--blacklist-bed", type=Path, default=None)

    parser.add_argument("--enrichr-top-genes", type=Path, default=None)
    parser.add_argument("--enrichr-ncg-genes", type=Path, default=None)

    parser.add_argument("--chunk-size", type=int, default=10000)
    parser.add_argument("--entropy-threshold", type=float, default=0.3)
    parser.add_argument("--min-recurrent-patients", type=int, default=10)
    parser.add_argument("--hotspot-merge-gap", type=int, default=200)

    parser.add_argument("--figure-group", default="all")

    parser.add_argument("--subgroup-min-patients", type=int, default=80)
    parser.add_argument("--subgroup-thresholds", default="0.40:0.90:0.02")
    parser.add_argument("--subgroup-target-main-groups", type=int, default=3)
    parser.add_argument("--subgroup-min-main-group-size", type=int, default=10)
    parser.add_argument("--subgroup-min-main-group-fraction", type=float, default=0.80)
    parser.add_argument("--subgroup-top-main-clusters", type=int, default=3)
    parser.add_argument("--subgroup-top-hotspots-per-group", type=int, default=5)
    parser.add_argument("--subgroup-signature-threshold", type=float, default=0.50)

    args = parser.parse_args()

    for folder in [
        args.data_external,
        args.data_interim,
        args.data_processed,
        args.results_tables,
        args.results_figures,
    ]:
        folder.mkdir(parents=True, exist_ok=True)

    commands = build_commands(args)

    for step in selected_steps(args):
        run_command(commands[step], dry_run=args.dry_run)

    print("\nPipeline finished.")


if __name__ == "__main__":
    main()
