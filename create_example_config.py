#!/usr/bin/env python3

"""
create_example_config.py

Create an example YAML configuration file for the funMotifs CRC ROSETTA
pipeline.

The current pipeline scripts are command-line driven, but this config file is
useful for documentation, reproducibility, and future wrapper development.

Example:
    python scripts/create_example_config.py --output config/config.example.yaml
"""

from __future__ import annotations

import argparse
from pathlib import Path


CONFIG_TEXT = """# Example configuration for the funMotifs CRC ROSETTA workflow.
# Edit these paths before running the pipeline on a new system.

project:
  name: funmotifs-crc-rosetta
  description: Identification of motif-disrupting variants in colorectal cancer using funMotifs and ROSETTA

paths:
  data_external: data/external
  data_interim: data/interim
  data_processed: data/processed
  results_tables: results/tables
  results_figures: results/figures
  results_logs: results/logs

inputs:
  # Colon motif annotation table used to prepare ROSETTA prediction input.
  colon_annotations: data/external/colon_annotations.tsv

  # Trained ROSETTA model objects.
  rosetta_rules: data/external/sig_rules_final_with_pretty_rules.rds
  rosetta_deployment: data/external/rosetta_deployment_info.rds

  # CRC variant BED file. The expected columns are defined in script 05.
  variant_bed: data/external/CRC-colon.section6.sorted.bed
  variant_bed_is_sorted: true

  # Motif position frequency matrix table.
  pfm: data/external/motifs_pfm.tsv

  # External interpretation files.
  gwas: data/external/gwas-association-CRC.tsv
  ncg: data/external/NCG_cancerdrivers_annotation_supporting_evidence.tsv
  gene_map: data/external/gene_id_to_symbol.tsv
  blacklist_bed: data/external/hg38-blacklist.v3.bed

  # Optional pathway enrichment exports, for example from Enrichr.
  enrichr_top_genes: data/external/KEGG_2026_table_T500_genes.txt
  enrichr_ncg_genes: data/external/KEGG_2026_T162_NCG.txt

parameters:
  # Script 02 and 03.
  chunk_size: 10000

  # Script 05.
  entropy_threshold: 0.30

  # Script 06.
  recurrent_min_patients: 10
  hotspot_merge_gap: 200

  # Script 08.
  figure_group: all

  # Script 09.
  subgroup_min_patients: 80
  subgroup_thresholds: "0.40:0.90:0.02"
  subgroup_target_main_groups: 3
  subgroup_min_main_group_size: 10
  subgroup_min_main_group_fraction: 0.80
  subgroup_top_main_clusters: 3
  subgroup_top_hotspots_per_group: 5
  subgroup_signature_threshold: 0.50

expected_outputs:
  rosetta_functional_table: data/processed/colon_rosetta_functional_only.tsv
  rosetta_functional_bed: data/processed/colon_rosetta_functional_motifs.sorted.bed
  motif_disruptions: results/tables/variant_overlap_disruption/motif_disrupting_variants.tsv
  recurrent_hotspots: results/tables/recurrent_hotspots_and_genes/recurrent_hotspot_regions.tsv
  gwas_besthit: results/tables/external_interpretation/crc_gwas_hotspot_besthit.tsv
  patient_subgroups: results/tables/patient_hotspot_subgroups/patient_jaccard_clusters.tsv
"""


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create example YAML config for the pipeline."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("config/config.example.yaml"),
        help="Output YAML path."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing config file."
    )

    args = parser.parse_args()

    if args.output.exists() and not args.overwrite:
        raise FileExistsError(
            f"File already exists: {args.output}. Use --overwrite to replace it."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(CONFIG_TEXT)

    print(f"Wrote example config: {args.output}")


if __name__ == "__main__":
    main()
