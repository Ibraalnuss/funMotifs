#!/usr/bin/env python3

"""
04_merge_rosetta_predictions.py

Merge per-chunk ROSETTA prediction outputs into final combined files.

After 03_predict_rosetta_colon.R writes one TSV and one BED file per chunk,
this script collects them and produces:

    colon_rosetta_functional_only.tsv      all functional motif rows
    colon_rosetta_predictions.tsv          full predictions (all classes), optional
    colon_rosetta_functional_motifs.bed    BED of functional motif coordinates
    colon_rosetta_functional_motifs.sorted.bed  same BED, sorted by coordinate
    rosetta_merge_summary.tsv              row counts and file paths

The BED sort is done in Python using pandas (no external unix sort needed).

Example:
    python scripts/04_merge_rosetta_predictions.py \\
      --input-dir results/tables/rosetta_chunk_predictions \\
      --out-functional data/processed/colon_rosetta_functional_only.tsv \\
      --out-full data/processed/colon_rosetta_predictions.tsv \\
      --out-bed data/processed/colon_rosetta_functional_motifs.bed \\
      --out-sorted-bed data/processed/colon_rosetta_functional_motifs.sorted.bed \\
      --out-summary results/tables/rosetta_merge_summary.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def count_data_rows(path: Path, has_header: bool = True) -> int:
    if not path.exists():
        return 0

    with path.open("r") as handle:
        n = sum(1 for _ in handle)

    if has_header and n > 0:
        return n - 1

    return n


def merge_tsv_files(files: list[Path], output: Path) -> int:
    """
    Merge TSV files with headers. Header is kept only once.
    """
    output.parent.mkdir(parents=True, exist_ok=True)

    if not files:
        print(f"No TSV files found for: {output}")
        return 0

    rows_written = 0

    with output.open("w") as out:
        for i, path in enumerate(files):
            with path.open("r") as handle:
                header = handle.readline()

                if i == 0:
                    out.write(header)

                for line in handle:
                    out.write(line)
                    rows_written += 1

    print(f"Saved merged TSV: {output}")
    print(f"Merged files: {len(files)}")
    print(f"Rows written: {rows_written:,}")

    return rows_written


def merge_bed_files(files: list[Path], output: Path) -> int:
    """
    Merge BED-like files without headers.
    """
    output.parent.mkdir(parents=True, exist_ok=True)

    if not files:
        print(f"No BED files found for: {output}")
        return 0

    rows_written = 0

    with output.open("w") as out:
        for path in files:
            with path.open("r") as handle:
                for line in handle:
                    if line.strip():
                        out.write(line)
                        rows_written += 1

    print(f"Saved merged BED: {output}")
    print(f"Merged BED files: {len(files)}")
    print(f"BED rows written: {rows_written:,}")

    return rows_written


def sort_bed(input_bed: Path, output_bed: Path) -> int:
    """
    Sort a BED-like file by chromosome, start, and end using pandas.

    This avoids requiring external unix sort/bedtools.
    """
    if not input_bed.exists() or input_bed.stat().st_size == 0:
        print(f"No BED rows to sort: {input_bed}")
        return 0

    df = pd.read_csv(input_bed, sep="\t", header=None)

    if df.empty:
        return 0

    if df.shape[1] < 3:
        raise ValueError(f"BED file must have at least 3 columns: {input_bed}")

    df["_chr_sort"] = df[0].astype(str).str.replace("chr", "", regex=False)

    chr_order = {
        **{str(i): i for i in range(1, 23)},
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 25,
    }

    df["_chr_sort"] = df["_chr_sort"].map(chr_order).fillna(999).astype(int)
    df["_start_sort"] = pd.to_numeric(df[1], errors="coerce").fillna(-1).astype(int)
    df["_end_sort"] = pd.to_numeric(df[2], errors="coerce").fillna(-1).astype(int)

    df = df.sort_values(["_chr_sort", "_start_sort", "_end_sort"])
    df = df.drop(columns=["_chr_sort", "_start_sort", "_end_sort"])

    output_bed.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_bed, sep="\t", header=False, index=False)

    print(f"Saved sorted BED: {output_bed}")
    print(f"Sorted BED rows: {len(df):,}")

    return len(df)


def find_files(input_dir: Path, pattern: str) -> list[Path]:
    return sorted(
        path for path in input_dir.glob(pattern)
        if path.is_file() and path.stat().st_size > 0
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge ROSETTA chunk predictions into final output files."
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Directory containing ROSETTA chunk prediction outputs."
    )
    parser.add_argument(
        "--out-functional",
        required=True,
        type=Path,
        help="Merged functional-only TSV output."
    )
    parser.add_argument(
        "--out-bed",
        required=True,
        type=Path,
        help="Merged functional motif BED output."
    )
    parser.add_argument(
        "--out-sorted-bed",
        required=True,
        type=Path,
        help="Sorted functional motif BED output."
    )
    parser.add_argument(
        "--out-summary",
        required=True,
        type=Path,
        help="Merge summary TSV output."
    )
    parser.add_argument(
        "--out-full",
        type=Path,
        default=None,
        help="Optional merged full prediction TSV output."
    )
    parser.add_argument(
        "--functional-pattern",
        default="*_functional_only.tsv",
        help="Glob pattern for functional-only chunk files."
    )
    parser.add_argument(
        "--full-pattern",
        default="*_predictions.tsv",
        help="Glob pattern for full prediction chunk files."
    )
    parser.add_argument(
        "--bed-pattern",
        default="*.bed",
        help="Glob pattern for BED files."
    )
    parser.add_argument(
        "--skip-full",
        action="store_true",
        help="Do not merge full prediction files."
    )

    args = parser.parse_args()

    if not args.input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {args.input_dir}")

    functional_files = find_files(args.input_dir, args.functional_pattern)
    full_files = find_files(args.input_dir, args.full_pattern)
    bed_files = find_files(args.input_dir, args.bed_pattern)

    # Avoid accidentally treating the final merged BED as an input if output is inside input_dir.
    bed_files = [
        path for path in bed_files
        if path.resolve() not in {
            args.out_bed.resolve(),
            args.out_sorted_bed.resolve(),
        }
    ]

    print("Input directory:", args.input_dir)
    print("Functional TSV files:", len(functional_files))
    print("Full prediction TSV files:", len(full_files))
    print("BED files:", len(bed_files))

    functional_rows = merge_tsv_files(functional_files, args.out_functional)

    full_rows = 0
    if not args.skip_full and args.out_full is not None:
        full_rows = merge_tsv_files(full_files, args.out_full)

    bed_rows = merge_bed_files(bed_files, args.out_bed)
    sorted_bed_rows = sort_bed(args.out_bed, args.out_sorted_bed)

    args.out_summary.parent.mkdir(parents=True, exist_ok=True)

    summary = pd.DataFrame([
        {
            "input_dir": str(args.input_dir),
            "functional_files": len(functional_files),
            "full_prediction_files": len(full_files),
            "bed_files": len(bed_files),
            "functional_rows": functional_rows,
            "full_prediction_rows": full_rows,
            "bed_rows": bed_rows,
            "sorted_bed_rows": sorted_bed_rows,
            "out_functional": str(args.out_functional),
            "out_full": str(args.out_full) if args.out_full else "",
            "out_bed": str(args.out_bed),
            "out_sorted_bed": str(args.out_sorted_bed),
        }
    ])

    summary.to_csv(args.out_summary, sep="\t", index=False)

    print(f"Saved merge summary: {args.out_summary}")


if __name__ == "__main__":
    main()
