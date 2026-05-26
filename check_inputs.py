#!/usr/bin/env python3

"""
check_inputs.py

Validate that the required input files and key columns are available before
running the funMotifs CRC ROSETTA pipeline.

This script does not run the analysis. It only checks file existence and
basic table structure.

Example:
    python scripts/check_inputs.py \
      --colon-annotations data/external/colon_annotations.tsv \
      --rules data/external/sig_rules_final_with_pretty_rules.rds \
      --deployment data/external/rosetta_deployment_info.rds \
      --variant-bed data/external/CRC-colon.section6.sorted.bed \
      --pfm data/external/motifs_pfm.tsv \
      --gwas data/external/gwas-association-CRC.tsv \
      --ncg data/external/NCG_cancerdrivers_annotation_supporting_evidence.tsv
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd


def ok(message: str) -> None:
    print(f"[OK] {message}")


def warn(message: str) -> None:
    print(f"[WARN] {message}")


def fail(message: str, errors: list[str]) -> None:
    print(f"[FAIL] {message}")
    errors.append(message)


def check_exists(path: Path | None, label: str, errors: list[str], required: bool = True) -> bool:
    if path is None:
        if required:
            fail(f"{label}: path not provided", errors)
            return False
        warn(f"{label}: path not provided, skipping optional check")
        return False

    if not path.exists():
        if required:
            fail(f"{label}: file not found: {path}", errors)
        else:
            warn(f"{label}: optional file not found: {path}")
        return False

    if path.is_file() and path.stat().st_size == 0:
        fail(f"{label}: file exists but is empty: {path}", errors)
        return False

    ok(f"{label}: {path}")
    return True


def read_header(path: Path, sep: str = "\t") -> list[str]:
    df = pd.read_csv(path, sep=sep, nrows=0)
    return list(df.columns)


def check_columns(
    path: Path,
    label: str,
    required_columns: list[str],
    errors: list[str],
    sep: str = "\t",
) -> bool:
    try:
        columns = read_header(path, sep=sep)
    except Exception as exc:
        fail(f"{label}: could not read header from {path}: {exc}", errors)
        return False

    missing = [col for col in required_columns if col not in columns]

    if missing:
        fail(f"{label}: missing columns {missing}", errors)
        print(f"       Available columns: {columns[:40]}{' ...' if len(columns) > 40 else ''}")
        return False

    ok(f"{label}: required columns present")
    return True


def check_bed_columns(path: Path, label: str, min_columns: int, errors: list[str]) -> bool:
    try:
        with path.open("r") as handle:
            for line in handle:
                if line.strip():
                    n_cols = len(line.rstrip("\n").split("\t"))
                    break
            else:
                fail(f"{label}: BED file has no data rows: {path}", errors)
                return False
    except Exception as exc:
        fail(f"{label}: could not read BED file {path}: {exc}", errors)
        return False

    if n_cols < min_columns:
        fail(f"{label}: expected at least {min_columns} columns, found {n_cols}", errors)
        return False

    ok(f"{label}: has at least {min_columns} BED columns")
    return True


def check_tool_available(tool: str, errors: list[str], required: bool = True) -> None:
    if shutil.which(tool):
        ok(f"command available: {tool}")
    else:
        if required:
            fail(f"required command not found in PATH: {tool}", errors)
        else:
            warn(f"optional command not found in PATH: {tool}")


def check_python_packages(errors: list[str]) -> None:
    packages = ["numpy", "pandas", "matplotlib", "scipy", "sklearn"]

    for pkg in packages:
        try:
            __import__(pkg)
            ok(f"Python package available: {pkg}")
        except Exception:
            fail(f"Python package missing: {pkg}", errors)


def check_r_files(args: argparse.Namespace, errors: list[str]) -> None:
    check_exists(args.rules, "ROSETTA rules RDS", errors, required=True)
    check_exists(args.deployment, "ROSETTA deployment RDS", errors, required=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate pipeline input files and key columns."
    )

    parser.add_argument("--colon-annotations", type=Path, default=None)
    parser.add_argument("--rules", type=Path, default=None)
    parser.add_argument("--deployment", type=Path, default=None)
    parser.add_argument("--variant-bed", type=Path, default=None)
    parser.add_argument("--functional-bed", type=Path, default=None)
    parser.add_argument("--pfm", type=Path, default=None)
    parser.add_argument("--gwas", type=Path, default=None)
    parser.add_argument("--ncg", type=Path, default=None)
    parser.add_argument("--gene-map", type=Path, default=None)
    parser.add_argument("--blacklist-bed", type=Path, default=None)
    parser.add_argument("--hotspot-table", type=Path, default=None)
    parser.add_argument("--event-table", type=Path, default=None)

    parser.add_argument(
        "--skip-tool-checks",
        action="store_true",
        help="Skip checking bedtools/Rscript availability."
    )
    parser.add_argument(
        "--skip-package-checks",
        action="store_true",
        help="Skip checking Python package imports."
    )

    args = parser.parse_args()
    errors: list[str] = []

    print("\nChecking command line tools")
    print("-" * 80)
    if not args.skip_tool_checks:
        check_tool_available("bedtools", errors, required=True)
        check_tool_available("Rscript", errors, required=True)
    else:
        warn("tool checks skipped")

    print("\nChecking Python packages")
    print("-" * 80)
    if not args.skip_package_checks:
        check_python_packages(errors)
    else:
        warn("Python package checks skipped")

    print("\nChecking core input files")
    print("-" * 80)

    if check_exists(args.colon_annotations, "colon annotation table", errors, required=True):
        check_columns(
            args.colon_annotations,
            "colon annotation table",
            ["mid", "chr", "motifstart", "motifend", "name", "strand"],
            errors,
        )

    check_r_files(args, errors)

    if check_exists(args.variant_bed, "CRC variant BED", errors, required=True):
        check_bed_columns(args.variant_bed, "CRC variant BED", 10, errors)

    if check_exists(args.pfm, "motif PFM table", errors, required=True):
        check_columns(args.pfm, "motif PFM table", ["name", "position", "allele", "freq"], errors)

    if check_exists(args.gwas, "CRC GWAS table", errors, required=True):
        check_columns(
            args.gwas,
            "CRC GWAS table",
            ["DISEASE/TRAIT", "CHR_ID", "CHR_POS", "MAPPED_GENE", "SNPS"],
            errors,
        )

    print("\nChecking optional inputs")
    print("-" * 80)

    if check_exists(args.ncg, "NCG annotation table", errors, required=False):
        check_columns(args.ncg, "NCG annotation table", ["symbol"], errors)

    if check_exists(args.gene_map, "gene ID to symbol map", errors, required=False):
        try:
            cols = read_header(args.gene_map)
            possible_gene = {"Gene", "gene_id", "ensembl_gene_id"}
            possible_symbol = {"Gene_Symbol", "gene_symbol", "hgnc_symbol", "external_gene_name"}
            if not possible_gene.intersection(cols) or not possible_symbol.intersection(cols):
                fail(
                    "gene map: expected one gene ID column and one symbol column",
                    errors,
                )
            else:
                ok("gene map: compatible columns present")
        except Exception as exc:
            fail(f"gene map: could not read header: {exc}", errors)

    if check_exists(args.blacklist_bed, "hg38 blacklist BED", errors, required=False):
        check_bed_columns(args.blacklist_bed, "hg38 blacklist BED", 3, errors)

    if check_exists(args.functional_bed, "functional motif BED", errors, required=False):
        check_bed_columns(args.functional_bed, "functional motif BED", 3, errors)

    if check_exists(args.hotspot_table, "recurrent hotspot table", errors, required=False):
        check_columns(
            args.hotspot_table,
            "recurrent hotspot table",
            ["Hotspot_ID", "Chr", "Hotspot_Start", "Hotspot_End", "Max_Unique_Patients"],
            errors,
        )

    if check_exists(args.event_table, "patient motif event table", errors, required=False):
        check_columns(
            args.event_table,
            "patient motif event table",
            ["Chr", "Motif_Start", "Motif_End", "Tumor_Sample_Barcode"],
            errors,
        )

    print("\nSummary")
    print("-" * 80)

    if errors:
        print(f"Input check failed with {len(errors)} problem(s):")
        for i, error in enumerate(errors, start=1):
            print(f"{i}. {error}")
        sys.exit(1)

    print("All required checks passed.")


if __name__ == "__main__":
    main()
