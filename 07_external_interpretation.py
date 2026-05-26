#!/usr/bin/env python3

"""
07_external_interpretation.py

Add external biological support to recurrent hotspot and gene results from
scripts 05 and 06. This is the annotation and validation step.

Steps:
  1. Filter GWAS Catalog rows to CRC-relevant traits (keyword matching).
  2. Compare recurrent hotspot regions to CRC GWAS loci; keep best hit per hotspot.
  3. Find exact position matches between somatic variants and GWAS SNPs.
  4. Optionally clean and annotate LD proxy hotspot overlap tables.
  5. Overlap top recurrent genes with NCG cancer driver annotations.
  6. Export gene lists for external pathway enrichment (e.g. Enrichr/KEGG).
  7. Optionally parse and plot Enrichr/KEGG pathway result tables as dotplots.

Main inputs:
    recurrent_hotspot_regions.tsv    from 06
    gene_prioritization_from_recurrent_loci.tsv    from 06
    motif_disrupting_variants.tsv    from 05
    gwas-association-CRC.tsv         GWAS Catalog download
    NCG_cancerdrivers_annotation_supporting_evidence.tsv    NCG database

Example:
    python scripts/07_external_interpretation.py \\
      --hotspots results/tables/recurrent_hotspots_and_genes/recurrent_hotspot_regions.tsv \\
      --genes results/tables/recurrent_hotspots_and_genes/gene_prioritization_from_recurrent_loci.tsv \\
      --disruptions results/tables/variant_overlap_disruption/motif_disrupting_variants.tsv \\
      --gwas data/external/gwas-association-CRC.tsv \\
      --outdir results/tables/external_interpretation
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


CRC_KEYWORDS = [
    "colorectal",
    "colon cancer",
    "rectal cancer",
    "colorectal cancer",
    "colorectal carcinoma",
    "colon carcinoma",
    "rectal carcinoma",
    "colorectal adenoma",
    "colorectal tumour",
    "colorectal tumor",
    "colorectal neoplasm",
    "colon neoplasm",
    "rectal neoplasm",
]


def clean_text(value) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def clean_chr(value) -> str:
    value = clean_text(value)
    if not value:
        return ""
    value = value.replace("CHR", "chr").replace("Chr", "chr")
    if not value.startswith("chr"):
        value = "chr" + value
    return value


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def save_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def trait_is_crc(trait: str) -> bool:
    trait = clean_text(trait).lower()
    return any(keyword in trait for keyword in CRC_KEYWORDS)


def split_gene_symbols(text: str) -> list[str]:
    text = clean_text(text)
    if not text:
        return []

    text = text.replace(" - ", ",")
    parts = []

    for block in text.split(","):
        for item in block.split(";"):
            item = item.strip()
            if item:
                parts.append(item)

    seen = set()
    out = []

    for gene in parts:
        if gene not in seen:
            seen.add(gene)
            out.append(gene)

    return out


def assign_evidence_tier(distance: int) -> str:
    if distance == 0:
        return "Tier 1: exact"
    if distance <= 1000:
        return "Tier 1: <=1kb"
    if distance <= 10000:
        return "Tier 2: <=10kb"
    if distance <= 50000:
        return "Tier 3: <=50kb"
    return ""


def trait_priority(trait: str) -> int:
    trait = clean_text(trait).lower()

    if "colorectal cancer" in trait or "colon cancer" in trait or "rectal cancer" in trait:
        return 1
    if "adenoma" in trait:
        return 2
    if "survival" in trait or "recurrence" in trait or "metastasis" in trait:
        return 3
    if "interaction" in trait or "pleiotropy" in trait:
        return 4
    return 5


def load_crc_gwas(gwas_file: Path, outdir: Path) -> pd.DataFrame:
    if not gwas_file.exists():
        raise FileNotFoundError(f"GWAS file not found: {gwas_file}")

    gwas = pd.read_csv(gwas_file, sep="\t", low_memory=False)

    required = [
        "DISEASE/TRAIT",
        "CHR_ID",
        "CHR_POS",
        "MAPPED_GENE",
        "SNPS",
        "PUBMEDID",
        "FIRST AUTHOR",
    ]
    missing = [col for col in required if col not in gwas.columns]
    if missing:
        raise ValueError(f"GWAS file is missing columns: {missing}")

    gwas = gwas.copy()
    gwas["DISEASE/TRAIT"] = gwas["DISEASE/TRAIT"].map(clean_text)
    gwas["CHR_ID"] = gwas["CHR_ID"].map(clean_text)
    gwas["CHR_POS"] = pd.to_numeric(gwas["CHR_POS"], errors="coerce")
    gwas["MAPPED_GENE"] = gwas["MAPPED_GENE"].map(clean_text)
    gwas["SNPS"] = gwas["SNPS"].map(clean_text)

    gwas = gwas.dropna(subset=["CHR_POS"]).copy()
    gwas = gwas[gwas["DISEASE/TRAIT"].map(trait_is_crc)].copy()

    if gwas.empty:
        raise RuntimeError("No CRC-relevant GWAS rows were found after filtering.")

    gwas["GWAS_Chr"] = gwas["CHR_ID"].map(clean_chr)
    gwas["GWAS_Pos"] = gwas["CHR_POS"].astype(int)

    keep_cols = [
        "DISEASE/TRAIT",
        "GWAS_Chr",
        "GWAS_Pos",
        "MAPPED_GENE",
        "SNPS",
        "PUBMEDID",
        "FIRST AUTHOR",
    ]
    optional_cols = [col for col in ["STUDY", "REGION"] if col in gwas.columns]

    out = gwas[keep_cols + optional_cols].drop_duplicates().reset_index(drop=True)

    out = out.rename(
        columns={
            "DISEASE/TRAIT": "GWAS_Trait",
            "MAPPED_GENE": "Mapped_Gene",
            "SNPS": "GWAS_SNP",
            "PUBMEDID": "PubMed_ID",
            "FIRST AUTHOR": "First_Author",
        }
    )

    if "STUDY" not in out.columns:
        out["STUDY"] = ""
    if "REGION" not in out.columns:
        out["REGION"] = ""

    save_tsv(out, outdir / "crc_gwas_filtered.tsv")
    return out


def compare_hotspots_to_gwas(
    hotspots: pd.DataFrame,
    gwas: pd.DataFrame,
    outdir: Path,
    max_distance: int,
    shortlist_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required_hotspot = ["Hotspot_ID", "Chr", "Hotspot_Start", "Hotspot_End", "Max_Unique_Patients"]
    missing = [col for col in required_hotspot if col not in hotspots.columns]
    if missing:
        raise ValueError(f"Hotspot table is missing columns: {missing}")

    rows = []

    for _, hs in hotspots.iterrows():
        chrom = clean_chr(hs["Chr"])
        start = int(hs["Hotspot_Start"])
        end = int(hs["Hotspot_End"])
        top_genes = set(split_gene_symbols(hs["Top_Genes"] if "Top_Genes" in hs.index else ""))

        same_chr = gwas[gwas["GWAS_Chr"] == chrom]

        for _, gw in same_chr.iterrows():
            pos = int(gw["GWAS_Pos"])

            if pos < start:
                distance = start - pos
            elif pos > end:
                distance = pos - end
            else:
                distance = 0

            if distance > max_distance:
                continue

            mapped_genes = split_gene_symbols(gw["Mapped_Gene"])
            concordant = any(gene in top_genes for gene in mapped_genes)

            rows.append(
                {
                    "Hotspot_ID": hs["Hotspot_ID"],
                    "Chr": chrom,
                    "Hotspot_Start": start,
                    "Hotspot_End": end,
                    "Hotspot_Width_bp": hs["Hotspot_Width_bp"] if "Hotspot_Width_bp" in hs.index else end - start + 1,
                    "Max_Unique_Patients": hs["Max_Unique_Patients"],
                    "Num_Loci_In_Hotspot": hs["Num_Loci_In_Hotspot"] if "Num_Loci_In_Hotspot" in hs.index else "",
                    "Top_TFs": hs["Top_TFs"] if "Top_TFs" in hs.index else "",
                    "Top_Genes": hs["Top_Genes"] if "Top_Genes" in hs.index else "",
                    "Representative_Locus": hs["Representative_Locus"] if "Representative_Locus" in hs.index else "",
                    "GWAS_Trait": gw["GWAS_Trait"],
                    "GWAS_SNP": gw["GWAS_SNP"],
                    "GWAS_Pos": gw["GWAS_Pos"],
                    "Mapped_Gene": gw["Mapped_Gene"],
                    "Distance_to_GWAS_bp": distance,
                    "Evidence_Tier": assign_evidence_tier(distance),
                    "Gene_Concordance": "yes" if concordant else "no",
                    "Trait_Priority": trait_priority(gw["GWAS_Trait"]),
                    "PubMed_ID": gw["PubMed_ID"],
                    "First_Author": gw["First_Author"],
                    "Study": gw["STUDY"],
                    "Region": gw["REGION"],
                }
            )

    full = pd.DataFrame(rows)

    if full.empty:
        save_tsv(full, outdir / "crc_gwas_hotspot_matches.tsv")
        save_tsv(full, outdir / "crc_gwas_hotspot_besthit.tsv")
        save_tsv(full, outdir / "crc_gwas_hotspot_besthit_shortlist.tsv")
        return full, full

    full = full.sort_values(
        ["Distance_to_GWAS_bp", "Trait_Priority", "Max_Unique_Patients"],
        ascending=[True, True, False],
    ).reset_index(drop=True)

    best = (
        full.sort_values(
            ["Distance_to_GWAS_bp", "Gene_Concordance", "Trait_Priority", "Max_Unique_Patients"],
            ascending=[True, False, True, False],
        )
        .drop_duplicates(subset=["Hotspot_ID"], keep="first")
        .sort_values(["Distance_to_GWAS_bp", "Max_Unique_Patients"], ascending=[True, False])
        .reset_index(drop=True)
    )

    shortlist = best.head(shortlist_n).copy()

    save_tsv(full, outdir / "crc_gwas_hotspot_matches.tsv")
    save_tsv(best, outdir / "crc_gwas_hotspot_besthit.tsv")
    save_tsv(shortlist, outdir / "crc_gwas_hotspot_besthit_shortlist.tsv")

    return full, best


def find_col(df: pd.DataFrame, candidates: list[str], required: bool = True) -> str | None:
    lower_map = {col.lower(): col for col in df.columns}

    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]

    if required:
        raise ValueError(f"Could not find any of these columns: {candidates}")

    return None


def find_exact_variant_gwas_matches(
    disruptions: pd.DataFrame,
    hotspots: pd.DataFrame,
    gwas: pd.DataFrame,
    outdir: Path,
    near_bp: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if disruptions.empty or hotspots.empty or gwas.empty:
        empty = pd.DataFrame()
        save_tsv(empty, outdir / "exact_variant_gwas_matches.tsv")
        save_tsv(empty, outdir / "exact_variant_gwas_matches_collapsed.tsv")
        save_tsv(empty, outdir / f"near_variant_gwas_matches_{near_bp}bp.tsv")
        return empty, empty

    variant_chr_col = find_col(disruptions, ["Variant_Chr", "Chr", "Chromosome", "chr"])
    variant_pos_col = find_col(
        disruptions,
        ["Variant_Start", "Start_Position", "Variant_Pos", "POS", "Position", "Motif_Start"],
    )

    sample_col = find_col(disruptions, ["Tumor_Sample_Barcode", "Sample", "Sample_ID"], required=False)
    ref_col = find_col(disruptions, ["Reference_Allele", "Ref", "ref"], required=False)
    alt_col = find_col(disruptions, ["Tumor_Seq_Allele2", "Alt", "alt", "Alternate_Allele"], required=False)
    gene_col = find_col(disruptions, ["Gene", "Nearest_Gene", "Target_Gene"], required=False)
    entropy_col = find_col(disruptions, ["Entropy", "Max_Entropy", "Disruption_Score"], required=False)
    motif_col = find_col(disruptions, ["Name", "Motif_Name", "Motif"], required=False)

    variants = pd.DataFrame(
        {
            "Chr": disruptions[variant_chr_col].map(clean_chr),
            "Variant_Pos": pd.to_numeric(disruptions[variant_pos_col], errors="coerce"),
            "Tumor_Sample_Barcode": disruptions[sample_col].astype(str) if sample_col else "",
            "Reference_Allele": disruptions[ref_col].astype(str) if ref_col else "",
            "Tumor_Seq_Allele2": disruptions[alt_col].astype(str) if alt_col else "",
            "Variant_Gene": disruptions[gene_col].astype(str) if gene_col else "",
            "Entropy": pd.to_numeric(disruptions[entropy_col], errors="coerce") if entropy_col else pd.NA,
            "Motif_Name": disruptions[motif_col].astype(str) if motif_col else "",
        }
    ).dropna(subset=["Chr", "Variant_Pos"])

    assigned = []

    for chrom in sorted(set(variants["Chr"]).intersection(set(hotspots["Chr"].map(clean_chr)))):
        vsub = variants[variants["Chr"] == chrom].copy()
        hsub = hotspots[hotspots["Chr"].map(clean_chr) == chrom].copy()

        for _, hs in hsub.iterrows():
            matched = vsub[
                (vsub["Variant_Pos"] >= int(hs["Hotspot_Start"]))
                & (vsub["Variant_Pos"] <= int(hs["Hotspot_End"]))
            ].copy()

            if matched.empty:
                continue

            for col in hotspots.columns:
                matched[col] = hs[col]

            assigned.append(matched)

    if not assigned:
        empty = pd.DataFrame()
        save_tsv(empty, outdir / "exact_variant_gwas_matches.tsv")
        save_tsv(empty, outdir / "exact_variant_gwas_matches_collapsed.tsv")
        save_tsv(empty, outdir / f"near_variant_gwas_matches_{near_bp}bp.tsv")
        return empty, empty

    vh = pd.concat(assigned, ignore_index=True)

    exact = vh.merge(
        gwas,
        left_on=["Chr", "Variant_Pos"],
        right_on=["GWAS_Chr", "GWAS_Pos"],
        how="inner",
    )

    near_rows = []

    for chrom in sorted(set(vh["Chr"]).intersection(set(gwas["GWAS_Chr"]))):
        vsub = vh[vh["Chr"] == chrom].copy()
        gsub = gwas[gwas["GWAS_Chr"] == chrom].copy()

        for _, gw in gsub.iterrows():
            pos = int(gw["GWAS_Pos"])
            nearby = vsub[(vsub["Variant_Pos"] >= pos - near_bp) & (vsub["Variant_Pos"] <= pos + near_bp)].copy()

            if nearby.empty:
                continue

            for col in gwas.columns:
                nearby[col] = gw[col]
            nearby["Distance_to_GWAS_bp"] = (nearby["Variant_Pos"] - pos).abs()
            near_rows.append(nearby)

    near = pd.concat(near_rows, ignore_index=True) if near_rows else pd.DataFrame()

    if not exact.empty:
        exact["Distance_to_GWAS_bp"] = 0

    save_tsv(exact, outdir / "exact_variant_gwas_matches.tsv")
    save_tsv(near, outdir / f"near_variant_gwas_matches_{near_bp}bp.tsv")

    if exact.empty:
        collapsed = pd.DataFrame()
    else:
        collapsed = (
            exact.groupby(
                [
                    "Chr",
                    "Variant_Pos",
                    "Reference_Allele",
                    "Tumor_Seq_Allele2",
                    "Hotspot_ID",
                    "Top_Genes",
                    "Top_TFs",
                    "GWAS_SNP",
                    "GWAS_Trait",
                    "Mapped_Gene",
                ],
                dropna=False,
            )
            .agg(
                Unique_Samples=("Tumor_Sample_Barcode", "nunique"),
                Samples=("Tumor_Sample_Barcode", lambda x: ",".join(sorted(set(map(str, x))))),
                Max_Unique_Patients=("Max_Unique_Patients", "max"),
                Hotspot_Width_bp=("Hotspot_Width_bp", "max"),
                Representative_Locus=("Representative_Locus", "first"),
                Max_Entropy=("Entropy", "max"),
            )
            .reset_index()
            .sort_values(["Unique_Samples", "Max_Unique_Patients"], ascending=[False, False])
        )

    save_tsv(collapsed, outdir / "exact_variant_gwas_matches_collapsed.tsv")
    return exact, collapsed


def clean_ld_overlap(ld_overlap_file: Path, hotspots: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    if ld_overlap_file is None or not ld_overlap_file.exists():
        return pd.DataFrame()

    raw = pd.read_csv(ld_overlap_file, sep="\t", header=None)

    raw.columns = [
        "Hotspot_Chr", "Hotspot_Start_bed", "Hotspot_End_bed", "Hotspot_ID", "Bed_Score", "Bed_Strand",
        "LD_Chr", "LD_Start_bed", "LD_End_bed", "LD_Name", "LD_Score", "LD_Strand",
        "Lead_GWAS_SNP", "Index_SNP", "Lead_GWAS_Coord", "Proxy_SNP", "Proxy_Coord",
        "R2", "Dprime", "Distance_to_Index", "Alleles", "MAF", "Correlated_Alleles",
        "FORGEdb", "RegulomeDB", "Functional_Class", "Source_File",
    ]

    raw["Proxy_Pos"] = raw["Proxy_Coord"].str.replace(r"chr[0-9XYM]+:", "", regex=True).astype(int)
    raw["Lead_GWAS_Pos"] = raw["Lead_GWAS_Coord"].str.replace(r"chr[0-9XYM]+:", "", regex=True).astype(int)

    merged = raw.merge(
        hotspots[
            [
                "Hotspot_ID",
                "Chr",
                "Hotspot_Start",
                "Hotspot_End",
                "Hotspot_Width_bp",
                "Num_Loci_In_Hotspot",
                "Max_Unique_Patients",
                "Top_TFs",
                "Top_Genes",
                "Representative_Locus",
            ]
        ],
        on="Hotspot_ID",
        how="left",
    )

    merged["GWAS_Overlap_Type"] = merged.apply(
        lambda row: "lead_SNP_direct" if row["Lead_GWAS_SNP"] == row["Proxy_SNP"] else "LD_proxy_overlap",
        axis=1,
    )

    final = merged[
        [
            "Hotspot_ID",
            "Chr",
            "Hotspot_Start",
            "Hotspot_End",
            "Hotspot_Width_bp",
            "Max_Unique_Patients",
            "Top_TFs",
            "Top_Genes",
            "Representative_Locus",
            "Lead_GWAS_SNP",
            "Lead_GWAS_Pos",
            "Proxy_SNP",
            "Proxy_Pos",
            "R2",
            "Dprime",
            "GWAS_Overlap_Type",
            "Distance_to_Index",
        ]
    ].copy()

    final = final.sort_values(["Lead_GWAS_SNP", "Hotspot_ID", "Proxy_Pos"]).reset_index(drop=True)
    save_tsv(final, outdir / "gwas_ld_overlap_hotspots_clean.tsv")
    return final


def overlap_top_genes_with_ncg(
    gene_table: pd.DataFrame,
    ncg_file: Path | None,
    outdir: Path,
    top_n: int,
) -> pd.DataFrame:
    if ncg_file is None or not ncg_file.exists() or gene_table.empty:
        empty = pd.DataFrame()
        save_tsv(empty, outdir / "top_genes_overlapping_NCG.tsv")
        save_text(outdir / "top_genes_overlapping_NCG_gene_symbols_only.txt", "")
        return empty

    ncg = pd.read_csv(ncg_file, sep="\t", low_memory=False)

    if "symbol" not in ncg.columns:
        raise ValueError("NCG file must contain column 'symbol'.")

    gene_col = None
    for candidate in ["Gene_Symbol", "Gene_Label", "Gene"]:
        if candidate in gene_table.columns:
            gene_col = candidate
            break

    if gene_col is None:
        raise ValueError("Gene table must contain Gene_Symbol, Gene_Label, or Gene.")

    rank_col = None
    for candidate in [
        "Unique_Patients",
        "Max_Locus_Unique_Patients",
        "Num_Total_Disrupting_Rows",
        "Unique_Recurrent_Loci",
    ]:
        if candidate in gene_table.columns:
            rank_col = candidate
            break

    genes = gene_table.copy()
    genes[gene_col] = genes[gene_col].map(clean_text)
    genes = genes[genes[gene_col] != ""].copy()

    if rank_col:
        genes = genes.sort_values(rank_col, ascending=False)

    genes = genes.drop_duplicates(subset=[gene_col]).head(top_n).copy()
    save_tsv(genes, outdir / f"top{top_n}_input_genes.tsv")

    ncg["symbol"] = ncg["symbol"].map(clean_text)
    ncg = ncg[ncg["symbol"] != ""].copy()

    keep_cols = [
        "symbol",
        "type",
        "organ_system",
        "primary_site",
        "cancer_type",
        "method",
        "coding_status",
        "cgc_annotation",
        "vogelstein_annotation",
        "saito_annotation",
        "NCG_oncogene",
        "NCG_tsg",
    ]
    keep_cols = [col for col in keep_cols if col in ncg.columns]
    ncg = ncg[keep_cols].drop_duplicates()

    overlap = genes.merge(ncg, left_on=gene_col, right_on="symbol", how="inner")

    if rank_col and rank_col in overlap.columns:
        overlap = overlap.sort_values(rank_col, ascending=False)

    save_tsv(overlap, outdir / f"top{top_n}_genes_overlapping_NCG.tsv")

    symbols = (
        overlap[gene_col]
        .dropna()
        .astype(str)
        .str.strip()
        .replace("", pd.NA)
        .dropna()
        .drop_duplicates()
        .tolist()
    )

    save_text(outdir / f"top{top_n}_genes_overlapping_NCG_gene_symbols_only.txt", "\n".join(symbols))
    return overlap


def export_pathway_inputs(gene_table: pd.DataFrame, ncg_overlap: pd.DataFrame, outdir: Path, top_n: int) -> None:
    gene_col = None

    for candidate in ["Gene_Symbol", "Gene_Label", "Gene"]:
        if candidate in gene_table.columns:
            gene_col = candidate
            break

    if gene_col is None or gene_table.empty:
        save_text(outdir / "pathway_input_top_genes.txt", "")
        save_text(outdir / "pathway_input_ncg_genes.txt", "")
        return

    genes = gene_table.copy()
    genes[gene_col] = genes[gene_col].map(clean_text)
    genes = genes[genes[gene_col] != ""].copy()

    rank_col = "Unique_Patients" if "Unique_Patients" in genes.columns else None
    if rank_col:
        genes = genes.sort_values(rank_col, ascending=False)

    top_genes = genes[gene_col].drop_duplicates().head(top_n).tolist()
    save_text(outdir / "pathway_input_top_genes.txt", "\n".join(top_genes))

    if not ncg_overlap.empty:
        ncg_gene_col = gene_col if gene_col in ncg_overlap.columns else "symbol"
        ncg_genes = (
            ncg_overlap[ncg_gene_col]
            .dropna()
            .astype(str)
            .str.strip()
            .replace("", pd.NA)
            .dropna()
            .drop_duplicates()
            .tolist()
        )
    else:
        ncg_genes = []

    save_text(outdir / "pathway_input_ncg_genes.txt", "\n".join(ncg_genes))


def read_table_auto(path: Path) -> pd.DataFrame:
    for sep in ["\t", ",", ";"]:
        try:
            df = pd.read_csv(path, sep=sep, low_memory=False)
            if df.shape[1] > 1:
                return df
        except Exception:
            continue
    raise ValueError(f"Could not parse table: {path}")


def split_gene_list(text: str) -> list[str]:
    text = clean_text(text)
    if not text:
        return []
    if ";" in text:
        parts = text.split(";")
    elif "," in text:
        parts = text.split(",")
    else:
        parts = [text]
    return [x.strip().upper() for x in parts if x.strip()]


def scale_sizes(values, min_size=120, max_size=1400):
    values = pd.Series(values).astype(float)

    if values.empty:
        return values

    vmin = values.min()
    vmax = values.max()

    if vmin == vmax:
        return pd.Series([0.5 * (min_size + max_size)] * len(values), index=values.index)

    return min_size + (values - vmin) * (max_size - min_size) / (vmax - vmin)


def format_enrichr_table(enrichr_file: Path, outdir: Path, label: str) -> pd.DataFrame:
    if enrichr_file is None or not enrichr_file.exists():
        return pd.DataFrame()

    df = read_table_auto(enrichr_file)

    required = ["Term", "Adjusted P-value", "Genes"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Enrichr file is missing columns: {missing}")

    df = df.copy()
    df["Adjusted P-value"] = pd.to_numeric(df["Adjusted P-value"], errors="coerce")
    df = df.dropna(subset=["Adjusted P-value"]).copy()

    rows = []

    for _, row in df.iterrows():
        genes = split_gene_list(row["Genes"])
        adjusted_p = float(row["Adjusted P-value"])
        rows.append(
            {
                "Term": clean_text(row["Term"]),
                "Adjusted_P_value": adjusted_p,
                "Minus_log10_Adjusted_P": -math.log10(max(adjusted_p, 1e-300)),
                "Genes": ";".join(genes),
                "Num_Genes": len(genes),
            }
        )

    out = pd.DataFrame(rows).sort_values(
        ["Minus_log10_Adjusted_P", "Num_Genes"],
        ascending=[False, False],
    ).reset_index(drop=True)

    save_tsv(out, outdir / f"{label}_pathway_table.tsv")
    return out


def make_pathway_dotplot(df: pd.DataFrame, out_fig: Path, title: str, top_n: int) -> None:
    if df.empty:
        return

    plot_df = df.head(top_n).copy().iloc[::-1].copy()
    plot_df["BubbleSize"] = scale_sizes(plot_df["Num_Genes"])

    fig_height = max(6, 0.55 * len(plot_df) + 2)

    fig, ax = plt.subplots(figsize=(12, fig_height))

    sc = ax.scatter(
        plot_df["Minus_log10_Adjusted_P"],
        range(len(plot_df)),
        s=plot_df["BubbleSize"],
        alpha=0.85,
    )

    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["Term"], fontsize=10)
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

    plt.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)

    hotspots = pd.read_csv(args.hotspots, sep="\t", low_memory=False)
    hotspots["Chr"] = hotspots["Chr"].map(clean_chr)

    gene_table = pd.read_csv(args.genes, sep="\t", low_memory=False) if args.genes and args.genes.exists() else pd.DataFrame()
    disruptions = pd.read_csv(args.disruptions, sep="\t", low_memory=False) if args.disruptions and args.disruptions.exists() else pd.DataFrame()

    gwas = load_crc_gwas(args.gwas, args.outdir)

    full_gwas, best_gwas = compare_hotspots_to_gwas(
        hotspots=hotspots,
        gwas=gwas,
        outdir=args.outdir,
        max_distance=args.gwas_window,
        shortlist_n=args.shortlist_n,
    )

    exact, exact_collapsed = find_exact_variant_gwas_matches(
        disruptions=disruptions,
        hotspots=hotspots,
        gwas=gwas,
        outdir=args.outdir,
        near_bp=args.near_bp,
    )

    ld_clean = clean_ld_overlap(
        ld_overlap_file=args.ld_overlap,
        hotspots=hotspots,
        outdir=args.outdir,
    )

    ncg_overlap = overlap_top_genes_with_ncg(
        gene_table=gene_table,
        ncg_file=args.ncg,
        outdir=args.outdir,
        top_n=args.top_n_genes,
    )

    export_pathway_inputs(
        gene_table=gene_table,
        ncg_overlap=ncg_overlap,
        outdir=args.outdir,
        top_n=args.top_n_genes,
    )

    pathway_tables = []

    if args.enrichr_top_genes:
        top_pathway = format_enrichr_table(args.enrichr_top_genes, args.outdir, "top_genes")
        pathway_tables.append(("top_genes", top_pathway))

        make_pathway_dotplot(
            top_pathway,
            args.outdir / "top_genes_pathway_dotplot.png",
            "Pathway enrichment for top recurrent genes",
            args.top_n_pathways,
        )

    if args.enrichr_ncg_genes:
        ncg_pathway = format_enrichr_table(args.enrichr_ncg_genes, args.outdir, "ncg_genes")
        pathway_tables.append(("ncg_genes", ncg_pathway))

        make_pathway_dotplot(
            ncg_pathway,
            args.outdir / "ncg_genes_pathway_dotplot.png",
            "Pathway enrichment for recurrent genes overlapping NCG",
            args.top_n_pathways,
        )

    summary_lines = [
        "External interpretation summary",
        "",
        f"CRC GWAS rows retained: {len(gwas):,}",
        f"Hotspot-GWAS matches within {args.gwas_window:,} bp: {len(full_gwas):,}",
        f"Best hotspot-GWAS rows: {len(best_gwas):,}",
        f"Exact variant-GWAS matches: {len(exact):,}",
        f"Collapsed exact variant-GWAS matches: {len(exact_collapsed):,}",
        f"Clean LD proxy overlap rows: {len(ld_clean):,}",
        f"NCG overlap rows: {len(ncg_overlap):,}",
    ]

    for label, table in pathway_tables:
        summary_lines.append(f"{label} pathway rows: {len(table):,}")

    save_text(args.outdir / "external_interpretation_summary.txt", "\n".join(summary_lines))

    print("\nDone.")
    print(f"Outputs saved in: {args.outdir}")
    print("\n".join(summary_lines))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run external interpretation for recurrent hotspots and genes."
    )

    parser.add_argument(
        "--hotspots",
        required=True,
        type=Path,
        help="Recurrent hotspot regions table from script 06."
    )
    parser.add_argument(
        "--genes",
        type=Path,
        default=None,
        help="Gene prioritization table from script 06."
    )
    parser.add_argument(
        "--disruptions",
        type=Path,
        default=None,
        help="Motif disrupting variants table from script 05."
    )
    parser.add_argument(
        "--gwas",
        required=True,
        type=Path,
        help="GWAS Catalog association table."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory."
    )
    parser.add_argument(
        "--ncg",
        type=Path,
        default=None,
        help="Optional NCG cancer driver annotation table."
    )
    parser.add_argument(
        "--ld-overlap",
        type=Path,
        default=None,
        help="Optional bedtools overlap table between hotspots and LDproxy BED."
    )
    parser.add_argument(
        "--enrichr-top-genes",
        type=Path,
        default=None,
        help="Optional Enrichr/KEGG export for top recurrent genes."
    )
    parser.add_argument(
        "--enrichr-ncg-genes",
        type=Path,
        default=None,
        help="Optional Enrichr/KEGG export for NCG overlapping recurrent genes."
    )
    parser.add_argument(
        "--gwas-window",
        type=int,
        default=50000,
        help="Maximum distance in bp between hotspot and GWAS SNP."
    )
    parser.add_argument(
        "--near-bp",
        type=int,
        default=100,
        help="Window for near variant-GWAS matches."
    )
    parser.add_argument(
        "--shortlist-n",
        type=int,
        default=20,
        help="Number of best GWAS-supported hotspots to keep in shortlist."
    )
    parser.add_argument(
        "--top-n-genes",
        type=int,
        default=500,
        help="Number of top recurrent genes for NCG/pathway input."
    )
    parser.add_argument(
        "--top-n-pathways",
        type=int,
        default=12,
        help="Number of pathway terms to show in dotplots."
    )

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
