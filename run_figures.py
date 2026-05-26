#!/usr/bin/env python3

"""
run_figures.py

Small convenience wrapper for 08_make_thesis_figures.py.

Use this when the tables are already generated and you only want to remake
figures.

Example:

    python scripts/run_figures.py --group all

    python scripts/run_figures.py --group gwas

    python scripts/run_figures.py \
      --recurrent-dir results/tables/recurrent_hotspots_and_genes \
      --external-dir results/tables/external_interpretation \
      --variant-dir results/tables/variant_overlap_disruption \
      --functional-table data/processed/colon_rosetta_functional_only.tsv \
      --outdir results/figures
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run thesis figure generation."
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
        help="Optional motif disrupting variants table."
    )
    parser.add_argument(
        "--functional-table",
        type=Path,
        default=Path("data/processed/colon_rosetta_functional_only.tsv"),
        help="Functional motif prediction table."
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

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print command without running it."
    )

    args = parser.parse_args()

    cmd = [
        sys.executable,
        str(SCRIPT_DIR / "08_make_thesis_figures.py"),
        "--group",
        args.group,
        "--recurrent-dir",
        str(args.recurrent_dir),
        "--external-dir",
        str(args.external_dir),
        "--variant-dir",
        str(args.variant_dir),
        "--outdir",
        str(args.outdir),
        "--dpi",
        str(args.dpi),
        "--total-motifs",
        str(args.total_motifs),
        "--functional-motifs",
        str(args.functional_motifs),
        "--top-loci",
        str(args.top_loci),
        "--top-labels",
        str(args.top_labels),
        "--top-gwas",
        str(args.top_gwas),
        "--top-gwas-schematic",
        str(args.top_gwas_schematic),
        "--top-tfs",
        str(args.top_tfs),
        "--top-genes",
        str(args.top_genes),
        "--top-pathways",
        str(args.top_pathways),
        "--min-patients-plot",
        str(args.min_patients_plot),
    ]

    if args.disruptions is not None:
        cmd.extend(["--disruptions", str(args.disruptions)])

    if args.functional_table is not None:
        cmd.extend(["--functional-table", str(args.functional_table)])

    if args.tf_background is not None:
        cmd.extend(["--tf-background", str(args.tf_background)])

    if args.require_gene_concordance:
        cmd.append("--require-gene-concordance")

    print("Running:")
    print(" ".join(cmd))

    if args.dry_run:
        return

    subprocess.run(cmd, check=True)

    print(f"\nFigure generation finished. Output directory: {args.outdir}")


if __name__ == "__main__":
    main()
