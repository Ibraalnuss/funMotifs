#!/usr/bin/env python3

"""
02_prepare_rosetta_colon_input.py

Convert a raw colon motif annotation table into the feature format expected by
the trained ROSETTA model, then optionally split the output into fixed-size
chunks for parallel R-based prediction.

What it does:
  1. Maps internal chromosome codes to standard names (e.g. 23 -> chrX).
  2. Converts cCRE, ChromHMM, replication timing domain, and signal columns
     into the binary/categorical predictor columns used during ROSETTA training.
  3. Applies Wang-style discretization to numothertfbinding and tfexpr.
  4. Optionally splits the processed table into chunks for 03_predict_rosetta_colon.R.

Required input columns:
    mid, chr, motifstart, motifend, name, strand
    (plus annotation columns: ccre, chromhmm, replidomain, dnase__seq, fantom,
     footprints, numothertfbinding, tfexpr)

Example:
    python scripts/02_prepare_rosetta_colon_input.py \\
      --input data/external/colon_annotations.tsv \\
      --output data/interim/colon_rosetta_input.tsv \\
      --chunk-dir data/interim/colon_rosetta_chunks \\
      --chunk-size 10000
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROSETTA_PREDICTORS = [
    "CTCF.only_CTCF.bound",
    "DNase.H3K4me3",
    "DNase.H3K4me3_CTCF.bound",
    "DNase.only",
    "Low.DNase",
    "PLS",
    "PLS_CTCF.bound",
    "dELS",
    "dELS_CTCF.bound",
    "pELS",
    "pELS_CTCF.bound",
    "Enh",
    "Quies",
    "Repr",
    "TSS",
    "dnase__seq",
    "fantom",
    "footprints",
    "DTZ",
    "ERD",
    "LRD",
    "UTZ",
    "numothertfbinding_wang",
    "tfexpr_wang",
]


CCRE_MAP = {
    "CTCF-only_CTCF-bound": "CTCF.only_CTCF.bound",
    "DNase-H3K4me3": "DNase.H3K4me3",
    "DNase-H3K4me3_CTCF-bound": "DNase.H3K4me3_CTCF.bound",
    "DNase-only": "DNase.only",
    "Low-DNase": "Low.DNase",
    "PLS": "PLS",
    "PLS_CTCF-bound": "PLS_CTCF.bound",
    "dELS": "dELS",
    "dELS_CTCF-bound": "dELS_CTCF.bound",
    "pELS": "pELS",
    "pELS_CTCF-bound": "pELS_CTCF.bound",
}


def safe_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def binarize_gt_zero(value):
    if pd.isna(value):
        return "0"
    return "1" if safe_float(value) > 0 else "0"


def safe_log10(value):
    if pd.isna(value):
        return 0.0

    value = safe_float(value)

    if value <= 0:
        return 0.0

    out = math.log10(value)

    if math.isnan(out) or math.isinf(out):
        return 0.0

    return out


def clip_numothertfbinding(value):
    if pd.isna(value):
        return 0.0
    return min(safe_float(value), 3.0)


def simplify_chromhmm_label(value):
    if pd.isna(value):
        return "NO"

    value = str(value)

    if value == "0":
        return "NO"
    if "Tss" in value or "Tx" in value or "BivFlnk" in value:
        return "TSS"
    if "Enh" in value:
        return "Enh"
    if "Repr" in value:
        return "Repr"
    if "Quies" in value or "Rpts" in value or "Het" in value:
        return "Quies"

    return "NO"


def map_chr_code(value):
    """
    Map internal chromosome codes to standard chromosome names.

    Expected:
    1..22 -> chr1..chr22
    23 -> chrX
    24 -> chrY
    25 -> chrM
    112 -> chr12
    """
    if pd.isna(value):
        raise ValueError("Missing chromosome code.")

    value = int(value)

    if 1 <= value <= 22:
        return f"chr{value}"
    if value == 23:
        return "chrX"
    if value == 24:
        return "chrY"
    if value == 25:
        return "chrM"
    if value == 112:
        return "chr12"

    raise ValueError(f"Unexpected chromosome code: {value}")


def wang_numothertfbinding(value, cut=2.55):
    return "1" if value <= cut else "2"


def wang_tfexpr(value, cut1=-0.1024527, cut2=0.1897982):
    if value <= cut1:
        return "1"
    if value <= cut2:
        return "2"
    return "3"


def make_empty_feature_frame(index):
    features = pd.DataFrame(index=index)

    for col in ROSETTA_PREDICTORS:
        features[col] = "0"

    return features


def preprocess_colon_annotations(
    raw_df,
    numother_cut=2.55,
    tfexpr_cut1=-0.1024527,
    tfexpr_cut2=0.1897982,
):
    features = make_empty_feature_frame(raw_df.index)

    if "ccre" in raw_df.columns:
        ccre = raw_df["ccre"].fillna("NO").astype(str)

        for raw_state, rosetta_col in CCRE_MAP.items():
            features[rosetta_col] = np.where(ccre == raw_state, "1", "0")

    if "chromhmm" in raw_df.columns:
        chrom = raw_df["chromhmm"].apply(simplify_chromhmm_label)

        for state in ["Enh", "Quies", "Repr", "TSS"]:
            features[state] = np.where(chrom == state, "1", "0")

    for col in ["dnase__seq", "fantom", "footprints"]:
        if col in raw_df.columns:
            features[col] = raw_df[col].apply(binarize_gt_zero)

    if "replidomain" in raw_df.columns:
        replidomain = raw_df["replidomain"].fillna("NO").astype(str)

        for state in ["DTZ", "ERD", "LRD", "UTZ"]:
            features[state] = np.where(replidomain == state, "1", "0")

    if "numothertfbinding" in raw_df.columns:
        numother = raw_df["numothertfbinding"].apply(clip_numothertfbinding)
        features["numothertfbinding_wang"] = numother.apply(
            lambda x: wang_numothertfbinding(x, numother_cut)
        )
    else:
        features["numothertfbinding_wang"] = "1"

    if "tfexpr" in raw_df.columns:
        tfexpr = raw_df["tfexpr"].apply(safe_log10)
        features["tfexpr_wang"] = tfexpr.apply(
            lambda x: wang_tfexpr(x, tfexpr_cut1, tfexpr_cut2)
        )
    else:
        features["tfexpr_wang"] = "2"

    return features[ROSETTA_PREDICTORS]


def build_rosetta_input(args):
    print(f"Reading colon annotations: {args.input}")

    raw_df = pd.read_csv(args.input, sep="\t", low_memory=False)
    print(f"Input shape: {raw_df.shape}")

    required = ["mid", "chr", "motifstart", "motifend", "name", "strand"]
    missing = [col for col in required if col not in raw_df.columns]

    if missing:
        raise ValueError("Missing required metadata columns: " + ", ".join(missing))

    meta = pd.DataFrame({
        "row_id": np.arange(1, len(raw_df) + 1),
        "mid": raw_df["mid"],
        "chr": raw_df["chr"].apply(map_chr_code),
        "motifstart": raw_df["motifstart"].astype(int),
        "motifend": raw_df["motifend"].astype(int),
        "name": raw_df["name"].astype(str),
        "strand": raw_df["strand"].astype(str),
    })

    features = preprocess_colon_annotations(
        raw_df,
        numother_cut=args.numother_cut,
        tfexpr_cut1=args.tfexpr_cut1,
        tfexpr_cut2=args.tfexpr_cut2,
    )

    out = pd.concat([meta, features], axis=1)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    print(f"Saved ROSETTA input: {args.output}")
    print(f"Output shape: {out.shape}")

    if args.chunk_dir:
        split_into_chunks(
            input_path=args.output,
            out_dir=args.chunk_dir,
            chunk_size=args.chunk_size,
            chunk_list=args.chunk_list,
        )


def split_into_chunks(input_path, out_dir, chunk_size, chunk_list=None):
    out_dir.mkdir(parents=True, exist_ok=True)

    chunk_paths = []
    rows_written = 0
    chunk_idx = 1

    with input_path.open("r") as handle:
        header = handle.readline()

        out_handle = None

        try:
            for line in handle:
                if rows_written % chunk_size == 0:
                    if out_handle is not None:
                        out_handle.close()

                    chunk_path = out_dir / f"colon_rosetta_input_chunk_{chunk_idx:03d}.tsv"
                    chunk_paths.append(chunk_path)

                    out_handle = chunk_path.open("w")
                    out_handle.write(header)

                    print(f"Writing chunk: {chunk_path}")
                    chunk_idx += 1

                out_handle.write(line)
                rows_written += 1

        finally:
            if out_handle is not None:
                out_handle.close()

    if chunk_list:
        chunk_list.parent.mkdir(parents=True, exist_ok=True)
        with chunk_list.open("w") as out:
            for path in chunk_paths:
                out.write(str(path) + "\n")
        print(f"Saved chunk list: {chunk_list}")

    print(f"Rows split: {rows_written:,}")
    print(f"Chunks written: {len(chunk_paths):,}")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare colon annotation data for ROSETTA prediction."
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Raw colon annotation table."
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output ROSETTA-ready TSV."
    )
    parser.add_argument(
        "--chunk-dir",
        type=Path,
        default=None,
        help="Optional output directory for split chunks."
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10000,
        help="Rows per chunk if --chunk-dir is provided."
    )
    parser.add_argument(
        "--chunk-list",
        type=Path,
        default=None,
        help="Optional text file listing all chunk paths."
    )
    parser.add_argument(
        "--numother-cut",
        type=float,
        default=2.55,
        help="Wang cutpoint for LR-preprocessed numothertfbinding."
    )
    parser.add_argument(
        "--tfexpr-cut1",
        type=float,
        default=-0.1024527,
        help="First Wang cutpoint for LR-preprocessed tfexpr."
    )
    parser.add_argument(
        "--tfexpr-cut2",
        type=float,
        default=0.1897982,
        help="Second Wang cutpoint for LR-preprocessed tfexpr."
    )

    args = parser.parse_args()
    build_rosetta_input(args)


if __name__ == "__main__":
    main()
