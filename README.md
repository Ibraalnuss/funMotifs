# funMotifs CRC — Identification of Motif-Disrupting Variants in Colorectal Cancer

This repository contains the analysis pipeline for the MSc thesis:

> **Identification of Motif-Disrupting Variants in Colorectal Cancer Using funMotifs**
> Ibraheem Al-Nussairi, Uppsala University, 2026

The pipeline predicts functionally active transcription factor (TF) motif instances in colon tissue, overlaps them with somatic variants from 1,063 CRC whole-genome samples, and prioritizes recurrently disrupted regulatory loci for interpretation.

---

## Table of Contents

- [Overview](#overview)
- [Repository structure](#repository-structure)
- [Dependencies](#dependencies)
- [Data requirements](#data-requirements)
- [Quick start](#quick-start)
- [Script reference](#script-reference)
  - [01 — Decision table preparation](#01--decision-table-preparation)
  - [02 — Prepare ROSETTA colon input](#02--prepare-rosetta-colon-input)
  - [03 — Predict with ROSETTA](#03--predict-with-rosetta)
  - [04 — Merge ROSETTA predictions](#04--merge-rosetta-predictions)
  - [05 — Variant overlap and disruption calling](#05--variant-overlap-and-disruption-calling)
  - [06 — Recurrent hotspots and gene prioritization](#06--recurrent-hotspots-and-gene-prioritization)
  - [07 — External interpretation](#07--external-interpretation)
  - [08 — Thesis figures](#08--thesis-figures)
  - [09 — Patient hotspot subgroups](#09--patient-hotspot-subgroups)
  - [run\_pipeline.py](#run_pipelinepy)
  - [run\_figures.py](#run_figurespy)
  - [check\_inputs.py](#check_inputspy)
  - [create\_example\_config.py](#create_example_configpy)
  - [summarize\_outputs.py](#summarize_outputspy)
- [Pipeline logic](#pipeline-logic)
- [Output files](#output-files)
- [Notes on reproducibility](#notes-on-reproducibility)

---

## Overview

The analysis proceeds in two main parts:

**Part 1 — Model building (scripts 01–04)**
A decision table of positive and negative motif instances is constructed and balanced, then used to train an interpretable ROSETTA rule-based model. The trained model is applied to colon tissue motif annotations to predict which motif instances are likely to be functionally active.

**Part 2 — Variant prioritization (scripts 05–09)**
Predicted functional motifs are overlapped with somatic CRC variants. SNPs are filtered by their impact on the motif position frequency matrix (PFM). The remaining candidate disruption events are aggregated into recurrent loci, hotspot regions, and gene/TF prioritization tables. External support from CRC GWAS loci, NCG cancer drivers, and pathway databases is added. Finally, patients are clustered into subgroups based on their shared hotspot disruption profiles.

---

## Repository structure

```
funmotifs_ibraheem/
│
├── scripts/                        # Numbered pipeline scripts (this repository)
│   ├── run_pipeline.py             # Wrapper: run any range of steps
│   ├── run_figures.py              # Wrapper: regenerate figures only
│   ├── check_inputs.py             # Validate required inputs and key columns
│   ├── create_example_config.py    # Write config/config.example.yaml
│   ├── summarize_outputs.py        # Summarize final pipeline outputs
│   ├── 01_prepare_decision_table.py
│   ├── 02_prepare_rosetta_colon_input.py
│   ├── 03_predict_rosetta_colon.R
│   ├── 04_merge_rosetta_predictions.py
│   ├── 05_variant_overlap_and_disruption.py
│   ├── 06_recurrent_hotspots_and_genes.py
│   ├── 07_external_interpretation.py
│   ├── 08_make_thesis_figures.py
│   └── 09_patient_hotspot_subgroups.py
│
├── data/
│   ├── external/                   # Raw input files (not tracked by git)
│   ├── interim/                    # Intermediate files generated during the run
│   └── processed/                  # Final prediction outputs
│
├── results/
│   ├── tables/                     # TSV output tables from each step
│   ├── figures/                    # PNG figures
│   └── logs/                       # Optional log files
│
├── config/                         # Optional configuration files
└── docs/                           # Notes and supplementary documentation
```

External data files and large intermediate outputs are not included in this repository. See [Data requirements](#data-requirements) for what you need to provide.

---

## Dependencies

### Python (≥ 3.9)

```
pandas
numpy
scikit-learn
matplotlib
scipy
pybedtools   (only needed for the intersect-mpra subcommand in script 01)
```

Install with conda (recommended):

```bash
conda create -n funmotifs python=3.10
conda activate funmotifs
conda install pandas numpy scikit-learn matplotlib scipy
conda install -c bioconda pybedtools bedtools
```

### R (≥ 4.2)

```r
install.packages("data.table")
# R.ROSETTA — install from GitHub if not on CRAN:
# devtools::install_github("komorowskilab/R.ROSETTA")
```

### System tools

`bedtools` must be available in your PATH for scripts 05 and 09.

---

## Data requirements

The following external files must be placed in `data/external/` before running the pipeline. They are not included in this repository.

| File | Source | Used by |
|---|---|---|
| `colon_annotations.tsv` | funMotifs colon tissue annotation output | 02 |
| `sig_rules_final_with_pretty_rules.rds` | Trained ROSETTA model (from thesis model training) | 03 |
| `rosetta_deployment_info.rds` | ROSETTA deployment metadata (pred_cols, factor_levels) | 03 |
| `wang_cutpoints_final.rds` | Wang discretization cutpoints (stored in deployment info) | 02, 03 |
| `CRC-colon.section6.sorted.bed` | CRC somatic variant BED (derived from MAF file) | 05 |
| `motifs_pfm.tsv` | Motif PFM table: columns name, position, allele, freq | 05 |
| `gwas-association-CRC.tsv` | GWAS Catalog download filtered to colorectal cancer | 07 |
| `NCG_cancerdrivers_annotation_supporting_evidence.tsv` | NCG cancer driver database | 07 |
| `KEGG_2026_table_T500_genes.txt` | Enrichr KEGG result for top recurrent genes (optional) | 07 |
| `KEGG_2026_T162_NCG.txt` | Enrichr KEGG result for NCG-overlapping genes (optional) | 07 |
| `hg38-blacklist.v3.bed` | ENCODE genomic blacklist (optional, used in script 09) | 09 |

The CRC variant BED file is derived from:

```
CRC-SW.Ensemble.1063_DNBSEQ.20210706.lite.maf.gz
```

This MAF file is not publicly available and is not included in the repository.

---

## Quick start

Create the example configuration file:

```bash
python scripts/create_example_config.py --output config/config.example.yaml
```

Before running the full pipeline, check that the required files and key columns are available:

```bash
python scripts/check_inputs.py \
  --colon-annotations data/external/colon_annotations.tsv \
  --rules data/external/sig_rules_final_with_pretty_rules.rds \
  --deployment data/external/rosetta_deployment_info.rds \
  --variant-bed data/external/CRC-colon.section6.sorted.bed \
  --pfm data/external/motifs_pfm.tsv \
  --gwas data/external/gwas-association-CRC.tsv \
  --ncg data/external/NCG_cancerdrivers_annotation_supporting_evidence.tsv
```

To run the full pipeline from ROSETTA prediction through patient subgroups:

```bash
python scripts/run_pipeline.py \
  --from rosetta_prepare --to subgroups \
  --colon-annotations data/external/colon_annotations.tsv \
  --rules data/external/sig_rules_final_with_pretty_rules.rds \
  --deployment data/external/rosetta_deployment_info.rds \
  --variant-bed data/external/CRC-colon.section6.sorted.bed \
  --pfm data/external/motifs_pfm.tsv \
  --gwas data/external/gwas-association-CRC.tsv \
  --ncg data/external/NCG_cancerdrivers_annotation_supporting_evidence.tsv
```

To check what would run without executing anything:

```bash
python scripts/run_pipeline.py --from rosetta_prepare --to subgroups --dry-run \
  --colon-annotations data/external/colon_annotations.tsv \
  ...
```

To regenerate figures only (after tables already exist):

```bash
python scripts/run_figures.py --group all
```

After a run, create a compact output summary for traceability:

```bash
python scripts/summarize_outputs.py \
  --data-processed data/processed \
  --results-tables results/tables \
  --out-prefix results/tables/pipeline_summary
```

---

## Script reference

### 01 — Decision table preparation

**File:** `scripts/01_prepare_decision_table.py`

A subcommand tool for preparing and evaluating the motif functionality decision table used for model training. Each subcommand does one task.

#### Subcommands

**`aggregate`** — Collapse cell-line-specific annotation columns (e.g. `hepg2___dnase__seq`, `k562___dnase__seq`) into one column per feature by taking the maximum across cells for numeric columns. Also binarizes signal columns (dnase, fantom, footprints), clips numothertfbinding, and log10-transforms tfexpr.

```bash
python scripts/01_prepare_decision_table.py aggregate \
  --input data/external/raw_annotations.tsv \
  --output data/interim/aggregated_features.tsv \
  --cells hepg2 k562
```

**`align`** — Align an aggregated feature table to a reference template, so that all tables have the same columns in the same order. Missing columns are filled with zeros. Optionally sets a label column to a fixed value (0 or 1).

```bash
python scripts/01_prepare_decision_table.py align \
  --input data/interim/aggregated_features.tsv \
  --template data/external/decision_table_template.tsv \
  --output data/interim/aligned_features.tsv \
  --label 1 \
  --label-column activity_score
```

**`balance`** — Subsample positives and negatives to equal counts (min of the two sizes), then shuffle rows. This avoids class imbalance during model training.

```bash
python scripts/01_prepare_decision_table.py balance \
  --positives data/interim/positives_aligned.tsv \
  --negatives data/interim/negatives_aligned.tsv \
  --output data/processed/decision_table_balanced.tsv
```

**`clean`** — Remove rows with obvious label conflicts. A labeled-positive row is removed if it has repressed chromatin state and no open chromatin or footprint support. A labeled-negative row is removed if it has a TSS state or active enhancer with open chromatin. This is optional but reduces label noise.

```bash
python scripts/01_prepare_decision_table.py clean \
  --input data/processed/decision_table_balanced.tsv \
  --output data/processed/decision_table_cleaned.tsv
```

**`evaluate`** — Split the table into train/test sets, fit a logistic regression model, and report accuracy, precision, recall, PR-AUC, and ROC-AUC for both splits. Also writes feature weights (beta coefficients and odds ratios) and a confusion matrix. Use `--permutation` to additionally compute permutation feature importance.

```bash
python scripts/01_prepare_decision_table.py evaluate \
  --table data/processed/decision_table_cleaned.tsv \
  --outdir results/tables/decision_table_evaluation \
  --name decision_table_v4
```

Key outputs in `--outdir`:

| File | Contents |
|---|---|
| `*_metrics.tsv` | Accuracy, precision, recall, PR-AUC, ROC-AUC for train and test |
| `*_confusion_matrix.tsv` | TN, FP, FN, TP on the test set |
| `*_predictions.csv` | y_true, y_proba, y_pred for every test row |
| `*_feature_weights.csv` | Beta coefficients and odds ratios per feature |
| `*_permutation_importance.csv` | Permutation importance scores (if --permutation) |

**`compare`** — Plot PR and ROC curves for two prediction CSV files side by side. Useful for comparing the baseline and improved decision tables.

```bash
python scripts/01_prepare_decision_table.py compare \
  --baseline results/tables/decision_table_evaluation/v1_predictions.csv \
  --new results/tables/decision_table_evaluation/v4_predictions.csv \
  --output results/figures/model_comparison.png \
  --baseline-label "Baseline (v1)" \
  --new-label "Improved (v4)"
```

**`intersect-mpra`** — Intersect a motif BED file with an MPRA region BED file to extract motifs that fall inside experimentally tested MPRA regions. Requires pybedtools.

```bash
python scripts/01_prepare_decision_table.py intersect-mpra \
  --motifs data/external/motifs.bed \
  --regions data/external/mpra_regions.bed \
  --output data/interim/motifs_in_mpra_regions.bed
```

---

### 02 — Prepare ROSETTA colon input

**File:** `scripts/02_prepare_rosetta_colon_input.py`

Converts a raw colon motif annotation table (output of funMotifs annotation) into the predictor format expected by the trained ROSETTA model, then optionally splits the result into fixed-size chunks for parallel prediction in R.

#### What it does

The raw colon annotation table contains columns like `ccre`, `chromhmm`, `replidomain`, `dnase__seq`, `numothertfbinding`, and `tfexpr`. This script converts those into the 24 binary/categorical predictor columns that the ROSETTA model was trained on:

- cCRE states (15 binary columns, e.g. `PLS`, `dELS`, `DNase.only`)
- ChromHMM simplified states (4 binary columns: `TSS`, `Enh`, `Repr`, `Quies`)
- Signal features binarized (`dnase__seq`, `fantom`, `footprints`)
- Replication timing domain states (4 binary columns: `DTZ`, `ERD`, `LRD`, `UTZ`)
- Wang-discretized continuous features (`numothertfbinding_wang`, `tfexpr_wang`)

Chromosome codes in the raw table (integers 1–25) are mapped to standard names (`chr1`–`chr22`, `chrX`, `chrY`, `chrM`).

#### Key functions

`preprocess_colon_annotations(raw_df)` — applies all feature conversions to a loaded DataFrame and returns a DataFrame with exactly the 24 ROSETTA predictor columns.

`build_rosetta_input(args)` — main function: reads the input, calls preprocessing, writes the output TSV, then optionally calls `split_into_chunks`.

`split_into_chunks(input_path, out_dir, chunk_size)` — splits the prepared TSV into files of `chunk_size` rows each, named `colon_rosetta_input_chunk_001.tsv`, `colon_rosetta_input_chunk_002.tsv`, etc.

#### Usage

```bash
python scripts/02_prepare_rosetta_colon_input.py \
  --input data/external/colon_annotations.tsv \
  --output data/interim/colon_rosetta_input.tsv \
  --chunk-dir data/interim/colon_rosetta_chunks \
  --chunk-size 10000 \
  --chunk-list data/interim/colon_rosetta_chunk_list.txt
```

#### Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--input` | yes | — | Raw colon annotation TSV |
| `--output` | yes | — | Output ROSETTA-ready TSV |
| `--chunk-dir` | no | — | If provided, split output into chunks here |
| `--chunk-size` | no | 10000 | Rows per chunk |
| `--chunk-list` | no | — | Text file listing all chunk paths (used by 03) |
| `--numother-cut` | no | 2.55 | Wang cutpoint for numothertfbinding |
| `--tfexpr-cut1` | no | -0.1024527 | First Wang cutpoint for tfexpr |
| `--tfexpr-cut2` | no | 0.1897982 | Second Wang cutpoint for tfexpr |

---

### 03 — Predict with ROSETTA

**File:** `scripts/03_predict_rosetta_colon.R`

Applies the trained ROSETTA rule model to every chunk produced by script 02. For each chunk, it coerces the predictor columns to the exact factor levels used during training, runs `predictClass()`, and writes the results.

#### What it does

For each chunk TSV:
1. Reads the chunk with `fread`.
2. Checks that all metadata and predictor columns are present.
3. Converts each predictor column to a factor with the levels stored in `rosetta_deployment_info.rds`.
4. Calls `predictClass(dt, rules, discrete=TRUE)` from R.ROSETTA.
5. Writes a full predictions TSV and a functional-only TSV for the chunk.
6. Appends functional motif rows to a combined BED file.

After all chunks are processed, a summary TSV with row counts and per-chunk timing is written.

#### ROSETTA deployment info

The `rosetta_deployment_info.rds` file must contain:
- `pred_cols`: character vector of predictor column names (must match what script 02 produces)
- `factor_levels`: named list where each element is the factor levels for one predictor column

The model itself is in `sig_rules_final_with_pretty_rules.rds` and was trained with R.ROSETTA using the settings: Johnson reducer, Modulo=TRUE, BRT=TRUE, BRTprec=0.99, Approximate=TRUE, Fraction=0.95.

#### Usage

```bash
Rscript scripts/03_predict_rosetta_colon.R \
  --chunk-dir data/interim/colon_rosetta_chunks \
  --rules data/external/sig_rules_final_with_pretty_ries.rds \
  --deployment data/external/rosetta_deployment_info.rds \
  --out-dir results/tables/rosetta_chunk_predictions
```

Alternatively, pass a chunk list file instead of a directory:

```bash
Rscript scripts/03_predict_rosetta_colon.R \
  --chunk-list data/interim/colon_rosetta_chunk_list.txt \
  --rules data/external/sig_rules_final_with_pretty_rules.rds \
  --deployment data/external/rosetta_deployment_info.rds \
  --out-dir results/tables/rosetta_chunk_predictions
```

#### Arguments

| Argument | Default | Description |
|---|---|---|
| `--chunk-dir` | — | Directory of chunk TSV files (alternative to --chunk-list) |
| `--chunk-list` | — | Text file listing chunk paths, one per line |
| `--rules` | `sig_rules_final_with_pretty_rules.rds` | Trained ROSETTA rules RDS |
| `--deployment` | `rosetta_deployment_info.rds` | Deployment metadata RDS |
| `--out-dir` | `rosetta_predictions` | Output directory |
| `--write-full` | TRUE | Write full prediction TSV per chunk |
| `--write-functional` | TRUE | Write functional-only TSV per chunk |
| `--write-bed` | TRUE | Append functional rows to combined BED |
| `--bed-name` | `colon_rosetta_functional_motifs.bed` | Name for the combined BED |
| `--summary-name` | `rosetta_prediction_summary.tsv` | Name for the per-chunk summary |

#### Outputs (per chunk)

```
rosetta_chunk_predictions/
├── colon_rosetta_input_chunk_001_predictions.tsv
├── colon_rosetta_input_chunk_001_functional_only.tsv
├── colon_rosetta_input_chunk_002_predictions.tsv
├── colon_rosetta_input_chunk_002_functional_only.tsv
├── ...
├── colon_rosetta_functional_motifs.bed
└── rosetta_prediction_summary.tsv
```

---

### 04 — Merge ROSETTA predictions

**File:** `scripts/04_merge_rosetta_predictions.py`

Collects all per-chunk output files from script 03 and merges them into final combined files. The BED merge also produces a coordinate-sorted version without requiring unix sort or bedtools.

#### Key functions

`merge_tsv_files(files, output)` — concatenates a list of TSV files with headers, writing the header only once.

`merge_bed_files(files, output)` — concatenates BED files (no header row).

`sort_bed(input_bed, output_bed)` — loads a BED into pandas, maps chromosome names to a numeric sort order, sorts by chromosome/start/end, and writes the sorted result.

`find_files(input_dir, pattern)` — returns all non-empty files in a directory matching a glob pattern.

#### Usage

```bash
python scripts/04_merge_rosetta_predictions.py \
  --input-dir results/tables/rosetta_chunk_predictions \
  --out-functional data/processed/colon_rosetta_functional_only.tsv \
  --out-full data/processed/colon_rosetta_predictions.tsv \
  --out-bed data/processed/colon_rosetta_functional_motifs.bed \
  --out-sorted-bed data/processed/colon_rosetta_functional_motifs.sorted.bed \
  --out-summary results/tables/rosetta_merge_summary.tsv
```

Add `--skip-full` to skip merging full prediction files (saves time and disk space if you only need functional motifs).

#### Outputs

| File | Contents |
|---|---|
| `colon_rosetta_functional_only.tsv` | All functional motif rows with metadata and prediction columns |
| `colon_rosetta_predictions.tsv` | All rows from all chunks including non-functional (optional) |
| `colon_rosetta_functional_motifs.bed` | Functional motifs in BED format: chr, start, end, name, strand, mid |
| `colon_rosetta_functional_motifs.sorted.bed` | Same, sorted by coordinate |
| `rosetta_merge_summary.tsv` | Row counts and file paths for all merged outputs |

---

### 05 — Variant overlap and disruption calling

**File:** `scripts/05_variant_overlap_and_disruption.py`

Overlaps the predicted functional colon motifs with CRC somatic variants and identifies candidates where a variant is likely to disrupt TF binding.

#### What it does

1. Sorts the functional motif BED by coordinate.
2. Optionally sorts the variant BED (skip with `--variant-bed-is-sorted`).
3. Runs `bedtools intersect -sorted -wo` to find all motif-variant overlaps.
4. Validates that the overlap file has the expected 17-column format.
5. For each overlapping SNP: looks up the reference and alternate allele frequencies at the variant position within the motif PFM. Retains the row if `ref_freq - alt_freq > entropy_threshold`. Non-SNP variants (indels, MNPs) are kept automatically.
6. Deduplicates, aggregates per patient, and writes summary tables.

#### Disruption logic

For an SNP at motif position `p`, the disruption score is:

```
disruption_score = PFM_freq(motif, p, ref_allele) - PFM_freq(motif, p, alt_allele)
```

A positive score means the variant replaces a high-frequency allele with a lower-frequency one — i.e. the mutation destabilizes TF binding at that position. The default threshold of 0.3 means only SNPs with a PFM frequency drop of at least 30% are retained.

Strand is taken into account: for minus-strand motifs, the position is computed as `motif_end - variant_start - 1`.

#### Key functions

`load_pfm_lookup(pfm_path)` — loads the PFM TSV and builds a dictionary keyed by `(motif_name, position, allele)` for fast frequency lookup.

`call_disrupting_variants(overlap_path, pfm_path, output_path, entropy_threshold)` — iterates over overlap rows, applies the disruption logic, and writes the filtered result.

`collapse_disrupted_loci(df_patient, outdir)` — groups by chr/start/end and counts unique patients per locus.

`cluster_motif_regions(df_regions, outdir, distance)` — groups motif regions within `distance` bp of each other into multi-motif clusters.

#### Usage

```bash
python scripts/05_variant_overlap_and_disruption.py \
  --functional-bed data/processed/colon_rosetta_functional_motifs.sorted.bed \
  --variant-bed data/external/CRC-colon.section6.sorted.bed \
  --pfm data/external/motifs_pfm.tsv \
  --outdir results/tables/variant_overlap_disruption \
  --variant-bed-is-sorted \
  --entropy-threshold 0.3
```

#### Outputs

```
results/tables/variant_overlap_disruption/
├── functional_motifs.sorted.bed
├── functional_motif_variant_overlaps.tsv      raw bedtools output
├── motif_disrupting_variants.tsv              filtered by disruption logic
└── cleaned/
    ├── motif_disrupting_variants_unique.tsv
    ├── unique_motif_region_per_patient.tsv     one row per motif x patient
    ├── unique_motif_regions.tsv
    ├── collapsed_disrupted_loci.tsv
    ├── top_tfs_by_unique_patients.tsv
    ├── top_motif_regions_by_unique_patients.tsv
    ├── multi_motif_clusters.tsv
    ├── analysis_summary.tsv
    └── overlap_disruption_funnel.tsv
```

---

### 06 — Recurrent hotspots and gene prioritization

**File:** `scripts/06_recurrent_hotspots_and_genes.py`

Counts unique patients per disrupted locus, filters by recurrence threshold, merges nearby loci into hotspot regions, and builds gene and TF prioritization tables.

#### What it does

1. Counts how many unique patients carry a disrupting variant at each chr/start/end locus.
2. Filters to loci seen in at least `--min-patients` patients.
3. Merges loci within `--merge-gap` bp of each other on the same chromosome into hotspot regions.
4. Links each hotspot to its associated genes and TFs (extracted from motif names using the `TF_MotifID` naming convention).
5. Builds an integrated locus/gene/TF table and separate gene-to-TF and TF-to-gene mapping tables.
6. Optionally annotates the top hotspots with functional feature summaries from the ROSETTA input chunks.
7. Optionally extracts non-TSS hotspot examples (loci where the disrupted motif is not in a TSS context).

#### Key functions

`motif_to_tf(motif_name)` — extracts the TF name from a motif ID by splitting on the last underscore (e.g. `ZNF384_MA1125.2` → `ZNF384`).

`build_recurrent_loci(disruptions_df, min_patients)` — counts patients per locus and returns filtered loci sorted by patient count.

`merge_loci_into_hotspots(recurrent_loci, merge_gap)` — greedy merge of nearby recurrent loci into hotspot intervals.

`build_gene_prioritization(disruptions_df, recurrent_loci)` — counts unique patients per gene linked to recurrent loci.

#### Usage

```bash
python scripts/06_recurrent_hotspots_and_genes.py \
  --disruptions results/tables/variant_overlap_disruption/motif_disrupting_variants.tsv \
  --outdir results/tables/recurrent_hotspots_and_genes \
  --min-patients 10 \
  --merge-gap 200 \
  --functional-bed data/processed/colon_rosetta_functional_motifs.sorted.bed \
  --functional-table data/processed/colon_rosetta_functional_only.tsv \
  --chunk-glob "data/interim/colon_rosetta_chunks/colon_rosetta_input_chunk_*.tsv"
```

#### Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--disruptions` | yes | — | motif_disrupting_variants.tsv from script 05 |
| `--outdir` | yes | — | Output directory |
| `--min-patients` | no | 10 | Minimum patients per locus to be called recurrent |
| `--merge-gap` | no | 200 | Max gap in bp for merging nearby loci into a hotspot |
| `--gene-map` | no | — | Optional gene ID to symbol mapping TSV |
| `--functional-bed` | no | — | Functional motif BED for feature annotation |
| `--functional-table` | no | — | Functional motif TSV for feature annotation |
| `--chunk-glob` | no | — | Glob pattern for ROSETTA input chunks |

#### Outputs

| File | Contents |
|---|---|
| `recurrent_loci_full.tsv` | All loci with patient counts (unfiltered) |
| `recurrent_loci_filtered.tsv` | Loci passing --min-patients threshold |
| `recurrent_hotspot_regions.tsv` | Merged hotspot regions with patient counts, TFs, genes |
| `recurrent_hotspot_regions.bed` | Same in BED format |
| `gene_prioritization_from_recurrent_loci.tsv` | Genes ranked by unique patient count |
| `pathway_input_genes_from_recurrent_loci.tsv` | Gene symbols for pathway enrichment |
| `integrated_recurrent_locus_table.tsv` | Full locus + gene + TF annotation table |
| `gene_to_top_tfs.tsv` | For each gene: the top associated TFs |
| `tf_to_top_genes.tsv` | For each TF: the top associated genes |
| `top_hotspot_functional_motif_summary.tsv` | Feature summaries for top hotspots (optional) |
| `non_tss_hotspot_examples.tsv` | Non-TSS hotspot examples (optional) |
| `recurrent_hotspots_and_genes_summary.txt` | Plain-text run summary |

---

### 07 — External interpretation

**File:** `scripts/07_external_interpretation.py`

Adds external biological support to the recurrent hotspot and gene results: CRC GWAS overlap, LD proxy matching, NCG cancer driver overlap, pathway enrichment input, and optional dotplots.

#### What it does

**GWAS overlap**
Filters the GWAS Catalog table to rows where the disease trait contains CRC-related keywords. For each recurrent hotspot, finds the nearest CRC GWAS locus and the one with the smallest p-value within a window (best hit). Also finds exact position matches where a somatic variant falls on the same base as a GWAS SNP.

**NCG overlap**
Takes the top N recurrent-locus genes (by unique patient count) and intersects them with the NCG cancer driver list. Returns the subset that have NCG support.

**Pathway enrichment**
Exports gene lists for external submission to Enrichr or KEGG. Optionally, if you have already run Enrichr and downloaded the result tables, this script parses them and generates bubble-style pathway dotplots.

#### Key functions

`trait_is_crc(trait)` — returns True if the trait string contains any CRC keyword (case-insensitive).

`filter_gwas_to_crc(gwas_df)` — applies `trait_is_crc` row-wise to subset the GWAS Catalog.

`compare_hotspots_to_gwas(hotspots_df, gwas_df, window)` — for each hotspot, computes distance to every CRC GWAS SNP on the same chromosome and keeps the best hit by p-value within `window` bp.

`find_exact_variant_gwas_matches(disruptions_df, gwas_df)` — finds rows where the somatic variant chromosome and position exactly match a CRC GWAS SNP position.

`overlap_genes_with_ncg(gene_df, ncg_df, top_n)` — takes the top N genes and returns those present in the NCG driver list.

`make_pathway_dotplot(pathway_df, output_path, title)` — generates a dotplot figure from an Enrichr/KEGG result table.

#### Usage

```bash
python scripts/07_external_interpretation.py \
  --hotspots results/tables/recurrent_hotspots_and_genes/recurrent_hotspot_regions.tsv \
  --genes results/tables/recurrent_hotspots_and_genes/gene_prioritization_from_recurrent_loci.tsv \
  --disruptions results/tables/variant_overlap_disruption/motif_disrupting_variants.tsv \
  --gwas data/external/gwas-association-CRC.tsv \
  --ncg data/external/NCG_cancerdrivers_annotation_supporting_evidence.tsv \
  --enrichr-top-genes data/external/KEGG_2026_table_T500_genes.txt \
  --enrichr-ncg-genes data/external/KEGG_2026_T162_NCG.txt \
  --outdir results/tables/external_interpretation
```

#### Outputs

| File | Contents |
|---|---|
| `crc_gwas_filtered.tsv` | GWAS Catalog rows matching CRC keywords |
| `crc_gwas_hotspot_besthit.tsv` | Best GWAS hit per hotspot (closest + smallest p-value) |
| `crc_gwas_hotspot_besthit_shortlist.tsv` | Hotspots with a GWAS hit within the window |
| `exact_variant_gwas_matches.tsv` | Somatic variants at exact GWAS SNP positions |
| `exact_variant_gwas_matches_collapsed.tsv` | Same, collapsed per locus |
| `near_variant_gwas_matches_100bp.tsv` | Variants within 100 bp of a GWAS SNP |
| `top500_genes_overlapping_NCG.tsv` | Top recurrent genes that are NCG cancer drivers |
| `pathway_input_top_genes.txt` | Gene list for Enrichr (top recurrent genes) |
| `pathway_input_ncg_genes.txt` | Gene list for Enrichr (NCG-overlapping genes) |
| `top_genes_pathway_table.tsv` | Parsed Enrichr/KEGG results for top genes |
| `ncg_genes_pathway_table.tsv` | Parsed Enrichr/KEGG results for NCG genes |
| `top_genes_pathway_dotplot.png` | Pathway dotplot for top genes |
| `ncg_genes_pathway_dotplot.png` | Pathway dotplot for NCG genes |
| `external_interpretation_summary.txt` | Plain-text run summary |

---

### 08 — Thesis figures

**File:** `scripts/08_make_thesis_figures.py`

Generates all thesis figures from the final analysis tables. All input paths are passed as arguments — nothing is hard-coded. Use `--group` to select one category or pass `all` to run everything.

#### Figure groups and outputs

| Group | Figures generated |
|---|---|
| `scoring` | `section5_scoring_summary.png`, `top_normalized_tfs_barplot.png` |
| `hotspots` | `manhattan_recurrent_loci.png`, `hotspot_landscape_manhattan_with_width.png` |
| `loci` | `top_integrated_recurrent_loci.png` |
| `gwas` | `gwas_hotspot_besthit_dotplot.png`, `crc_gwas_hotspot_direct_overlap_schematic.png` |
| `tf_gene` | `top_tfs_by_unique_patients_detailed.png`, `tf_gene_intersection_bubble_plot.png` |
| `pathways` | `top_genes_pathway_dotplot.png`, `ncg_genes_pathway_dotplot.png` |

Manhattan-style plots use hg38 chromosome lengths and plot recurrent loci by genomic position with patient count on the y-axis. The combined functional BED is used for scoring summaries.

#### Usage

```bash
python scripts/08_make_thesis_figures.py \
  --group all \
  --recurrent-dir results/tables/recurrent_hotspots_and_genes \
  --external-dir results/tables/external_interpretation \
  --variant-dir results/tables/variant_overlap_disruption \
  --functional-table data/processed/colon_rosetta_functional_only.tsv \
  --outdir results/figures \
  --dpi 300
```

Individual groups:

```bash
python scripts/08_make_thesis_figures.py --group gwas \
  --recurrent-dir results/tables/recurrent_hotspots_and_genes \
  --external-dir results/tables/external_interpretation \
  --outdir results/figures
```

#### Key arguments

| Argument | Default | Description |
|---|---|---|
| `--group` | `all` | Figure group(s) to generate |
| `--recurrent-dir` | — | Output directory from script 06 |
| `--external-dir` | — | Output directory from script 07 |
| `--variant-dir` | — | Output directory from script 05 |
| `--functional-table` | — | colon_rosetta_functional_only.tsv from script 04 |
| `--outdir` | `results/figures` | Where to save figures |
| `--dpi` | 300 | Figure resolution |
| `--top-loci` | 15 | Loci to show in barplots |
| `--top-tfs` | 20 | TFs to show in TF plots |
| `--top-genes` | 12 | Genes to show in gene plots |
| `--top-pathways` | 12 | Pathways to show in dotplots |
| `--total-motifs` | 0 | Total motif count for scoring summary annotation |
| `--functional-motifs` | 0 | Functional motif count for scoring summary annotation |

---

### 09 — Patient hotspot subgroups

**File:** `scripts/09_patient_hotspot_subgroups.py`

Groups CRC patients into subgroups based on their shared recurrent hotspot disruption profiles using Jaccard distance and hierarchical clustering.

#### What it does

1. Optionally removes hotspots overlapping a genomic blacklist BED.
2. Filters hotspots to those seen in at least `--min-patients` patients (default: 80).
3. Uses `bedtools intersect` to find which patients carry events at each selected hotspot.
4. Builds a binary patient x hotspot matrix (1 = patient has a disrupting event at that hotspot).
5. Computes pairwise Jaccard distances and runs average-linkage hierarchical clustering.
6. Scans clustering thresholds in steps of 0.02 across a user-defined range. For each threshold, counts how many clusters have at least `--min-main-group-size` patients. Chooses the threshold that gives close to `--target-main-groups` interpretable groups.
7. Assigns patients to clusters, relabels them C1, C2, ... by descending size.
8. Summarizes hotspot frequencies and top genes/TFs per cluster.
9. Exports signature and full gene lists per cluster for Enrichr.
10. Generates five figures (see below).

#### Key functions

`build_patient_hotspot_matrix(intersections_df)` — pivots the bedtools intersect output into a binary patient x hotspot DataFrame.

`jaccard_clustering(matrix_df)` — computes pairwise Jaccard distances using scipy.spatial.distance.pdist with metric='jaccard', then runs average-linkage hierarchical clustering.

`scan_thresholds(linkage_matrix, thresholds, min_group_size, min_fraction)` — evaluates each threshold in the scan range and returns the threshold scan table.

`assign_clusters(linkage_matrix, threshold)` — cuts the dendrogram at a given threshold and relabels clusters C1, C2, ... by size.

`summarize_cluster_profiles(patient_clusters, hotspot_matrix, hotspot_table)` — computes hotspot disruption frequency per cluster.

#### Usage

```bash
python scripts/09_patient_hotspot_subgroups.py \
  --hotspot-table results/tables/recurrent_hotspots_and_genes/recurrent_hotspot_regions.tsv \
  --event-table results/tables/variant_overlap_disruption/cleaned/unique_motif_region_per_patient.tsv \
  --outdir results/tables/patient_hotspot_subgroups \
  --figdir results/figures/patient_hotspot_subgroups \
  --min-patients 80 \
  --thresholds 0.40:0.90:0.02 \
  --target-main-groups 3
```

#### Key arguments

| Argument | Default | Description |
|---|---|---|
| `--hotspot-table` | — | recurrent_hotspot_regions.tsv from script 06 |
| `--event-table` | — | unique_motif_region_per_patient.tsv from script 05 |
| `--outdir` | — | Output directory for tables |
| `--figdir` | — | Output directory for figures |
| `--min-patients` | 80 | Minimum patients for a hotspot to be included |
| `--thresholds` | `0.40:0.90:0.02` | Threshold scan range: start:stop:step |
| `--target-main-groups` | 3 | Desired number of main patient groups |
| `--min-main-group-size` | 10 | Minimum patients for a group to count as "main" |
| `--min-main-group-fraction` | 0.80 | Fraction of all patients that must be in main groups |
| `--top-hotspots-per-group` | 5 | Hotspots to show per group in profile plots |
| `--signature-threshold` | 0.50 | Minimum disruption frequency in a cluster for "signature" hotspots |
| `--blacklist-bed` | — | Optional blacklist BED for filtering hotspots |

#### Outputs (tables)

| File | Contents |
|---|---|
| `selected_hotspots_ge80.tsv` | Hotspots passing --min-patients filter |
| `patient_by_hotspot_binary_matrix.tsv` | Patient x hotspot binary presence matrix |
| `jaccard_threshold_scan.tsv` | Main group counts at each threshold |
| `patient_jaccard_clusters.tsv` | Cluster assignment per patient |
| `patient_cluster_sizes.tsv` | Patient count per cluster |
| `cluster_hotspot_profiles.tsv` | Hotspot disruption frequency per cluster |
| `cluster_gene_tf_summary.tsv` | Top genes and TFs per cluster |
| `gene_lists_for_pathway_enrichment/` | Signature and full gene lists per cluster |
| `patient_hotspot_subgroups_summary.txt` | Plain-text run summary |

#### Outputs (figures)

| File | Description |
|---|---|
| `jaccard_threshold_scan.png` | Line plot: number of main groups across thresholds |
| `patient_group_sizes.png` | Barplot: patients per cluster |
| `cluster_hotspot_frequency_heatmap.png` | Heatmap: hotspot disruption frequency per cluster |
| `patient_jaccard_dendrogram.png` | Hierarchical clustering dendrogram |
| `clear_subgroup_hotspot_profiles.png` | Profile barplots for top hotspots per cluster |

---

### run\_pipeline.py

**File:** `scripts/run_pipeline.py`

A convenience wrapper that calls the numbered pipeline scripts in the correct order. It does not contain analysis logic — it just builds and runs the right commands.

**Step names and their corresponding scripts:**

| Step name | Script |
|---|---|
| `decision_table` | 01_prepare_decision_table.py (evaluate subcommand) |
| `rosetta_prepare` | 02_prepare_rosetta_colon_input.py |
| `rosetta_predict` | 03_predict_rosetta_colon.R |
| `rosetta_merge` | 04_merge_rosetta_predictions.py |
| `variant_overlap` | 05_variant_overlap_and_disruption.py |
| `recurrent` | 06_recurrent_hotspots_and_genes.py |
| `external` | 07_external_interpretation.py |
| `figures` | 08_make_thesis_figures.py |
| `subgroups` | 09_patient_hotspot_subgroups.py |

```bash
# Run a single step
python scripts/run_pipeline.py --step variant_overlap --variant-bed ... --pfm ...

# Run from step A to step B
python scripts/run_pipeline.py --from rosetta_prepare --to recurrent \
  --colon-annotations ... --rules ... --deployment ... --variant-bed ... --pfm ...

# Check commands without running
python scripts/run_pipeline.py --from rosetta_prepare --to subgroups --dry-run \
  --colon-annotations ... --rules ... --deployment ... --variant-bed ... --pfm ...
```

---

### run\_figures.py

**File:** `scripts/run_figures.py`

A small wrapper around script 08 for when the analysis tables already exist and you only want to regenerate or update figures.

```bash
# Regenerate all figures with defaults
python scripts/run_figures.py --group all

# Regenerate only GWAS figures
python scripts/run_figures.py --group gwas

# Use non-default input directories
python scripts/run_figures.py \
  --group all \
  --recurrent-dir results/tables/recurrent_hotspots_and_genes \
  --external-dir results/tables/external_interpretation \
  --variant-dir results/tables/variant_overlap_disruption \
  --outdir results/figures
```


---

### check\_inputs.py

**File:** `scripts/check_inputs.py`

Validates that the main pipeline inputs are present and that the most important table columns are available before running the workflow. This script does not run any analysis. It is intended as a fast pre-flight check so that missing files, wrong paths, or incompatible column names are caught early.

#### What it checks

- Required command-line tools: `bedtools` and `Rscript`
- Required Python packages: `numpy`, `pandas`, `matplotlib`, `scipy`, and `sklearn`
- Core input files: colon annotation table, ROSETTA rules, ROSETTA deployment metadata, CRC variant BED, motif PFM table, and GWAS table
- Optional inputs when supplied: NCG table, gene map, blacklist BED, functional motif BED, recurrent hotspot table, and patient motif event table

#### Required column checks

| Input | Required columns or structure |
|---|---|
| Colon annotations | `mid`, `chr`, `motifstart`, `motifend`, `name`, `strand` |
| CRC variant BED | At least 10 BED-like columns |
| Motif PFM table | `name`, `position`, `allele`, `freq` |
| GWAS table | `DISEASE/TRAIT`, `CHR_ID`, `CHR_POS`, `MAPPED_GENE`, `SNPS` |
| NCG table | `symbol` |
| Hotspot table | `Hotspot_ID`, `Chr`, `Hotspot_Start`, `Hotspot_End`, `Max_Unique_Patients` |
| Patient event table | `Chr`, `Motif_Start`, `Motif_End`, `Tumor_Sample_Barcode` |

#### Usage

```bash
python scripts/check_inputs.py \
  --colon-annotations data/external/colon_annotations.tsv \
  --rules data/external/sig_rules_final_with_pretty_rules.rds \
  --deployment data/external/rosetta_deployment_info.rds \
  --variant-bed data/external/CRC-colon.section6.sorted.bed \
  --pfm data/external/motifs_pfm.tsv \
  --gwas data/external/gwas-association-CRC.tsv \
  --ncg data/external/NCG_cancerdrivers_annotation_supporting_evidence.tsv
```

Optional checks can be added when those files are available:

```bash
python scripts/check_inputs.py \
  --colon-annotations data/external/colon_annotations.tsv \
  --rules data/external/sig_rules_final_with_pretty_rules.rds \
  --deployment data/external/rosetta_deployment_info.rds \
  --variant-bed data/external/CRC-colon.section6.sorted.bed \
  --pfm data/external/motifs_pfm.tsv \
  --gwas data/external/gwas-association-CRC.tsv \
  --ncg data/external/NCG_cancerdrivers_annotation_supporting_evidence.tsv \
  --gene-map data/external/gene_id_to_symbol.tsv \
  --blacklist-bed data/external/hg38-blacklist.v3.bed
```

Use `--skip-tool-checks` or `--skip-package-checks` only if you want to validate file structure without checking the local software environment.

---

### create\_example\_config.py

**File:** `scripts/create_example_config.py`

Creates a documented example YAML configuration file at `config/config.example.yaml`. The current workflow is still command-line driven, but this file provides a clear overview of expected paths, parameters, and major outputs. It is useful for documentation, reproducibility, and future development if the pipeline is later changed to read from a config file directly.

#### Usage

```bash
python scripts/create_example_config.py --output config/config.example.yaml
```

If the file already exists and you want to replace it:

```bash
python scripts/create_example_config.py \
  --output config/config.example.yaml \
  --overwrite
```

#### Contents of the generated config

The generated YAML contains four main sections:

| Section | Purpose |
|---|---|
| `project` | Project name and short description |
| `paths` | Standard locations for external data, interim files, processed data, tables, figures, and logs |
| `inputs` | Expected input files, including colon annotations, ROSETTA objects, CRC variant BED, PFM table, GWAS table, NCG table, and optional enrichment files |
| `parameters` | Key pipeline parameters such as chunk size, entropy threshold, recurrence threshold, hotspot merge gap, and subgroup clustering settings |
| `expected_outputs` | Main output files that should exist after a successful run |

---

### summarize\_outputs.py

**File:** `scripts/summarize_outputs.py`

Collects row counts and key values from the main pipeline outputs and writes a compact summary table and text report. This is mainly for thesis traceability and final quality control. It does not create new biological results; it summarizes what the pipeline already produced.

#### Usage

```bash
python scripts/summarize_outputs.py \
  --data-processed data/processed \
  --results-tables results/tables \
  --out-prefix results/tables/pipeline_summary
```

#### Outputs

| File | Contents |
|---|---|
| `pipeline_summary.tsv` | Machine-readable summary with section, metric, value, and source file |
| `pipeline_summary.txt` | Human-readable report grouped by analysis section |

#### What it summarizes

The summary includes values such as:

- Number of predicted functional motifs
- Number of motif-variant overlaps
- Number of retained motif-disrupting events
- Number of unique patients with retained disruptions
- Number of recurrent loci and hotspot regions
- Number of prioritized genes
- Number of GWAS-supported hotspots
- Number of exact somatic variant and GWAS matches
- Number of NCG-overlapping genes
- Number of patient subgroups and the chosen Jaccard threshold

---

## Pipeline logic

The overall analysis flow is:

```
create_example_config.py  → optional example config
check_inputs.py           → pre-flight validation
        │
        ▼
colon_annotations.tsv
        │
        ▼
02  preprocess → predictor columns, chunk split
        │
        ▼
03  ROSETTA predict (R) → per-chunk functional predictions
        │
        ▼
04  merge chunks → colon_rosetta_functional_motifs.sorted.bed
        │
        ▼
05  bedtools intersect + PFM filter → motif_disrupting_variants.tsv
        │
        ▼
06  recurrence aggregation → recurrent_hotspot_regions.tsv
                                 gene_prioritization_from_recurrent_loci.tsv
        │
        ▼
07  external annotation → GWAS overlap, NCG overlap, pathway enrichment
        │
        ▼
08  figures → thesis figures
        │
        ▼
09  patient subgroups → Jaccard clustering, subgroup figures
        │
        ▼
summarize_outputs.py → compact run summary
```

The decision table evaluation (script 01) is independent of this flow and documents the model training process.

---

## Output files

The most important downstream-facing output files are:

| File | Location | Description |
|---|---|---|
| `colon_rosetta_functional_motifs.sorted.bed` | `data/processed/` | All functional motifs in colon tissue |
| `motif_disrupting_variants.tsv` | `results/tables/variant_overlap_disruption/` | All retained candidate disruption events |
| `unique_motif_region_per_patient.tsv` | `.../cleaned/` | One row per patient x motif region |
| `recurrent_hotspot_regions.tsv` | `results/tables/recurrent_hotspots_and_genes/` | Prioritized recurrent hotspot regions |
| `gene_prioritization_from_recurrent_loci.tsv` | same | Genes ranked by patient count |
| `crc_gwas_hotspot_besthit.tsv` | `results/tables/external_interpretation/` | GWAS-supported hotspots |
| `exact_variant_gwas_matches.tsv` | same | Variants landing on GWAS SNP positions |
| `top500_genes_overlapping_NCG.tsv` | same | Recurrent genes with NCG cancer driver support |
| `patient_jaccard_clusters.tsv` | `results/tables/patient_hotspot_subgroups/` | Patient cluster assignments |
| `pipeline_summary.tsv` | `results/tables/` | Compact summary of the main pipeline outputs |
| `pipeline_summary.txt` | `results/tables/` | Human-readable summary report for traceability |

---

## Notes on reproducibility

The numbered scripts in `scripts/` represent the final, cleaned version of the analysis workflow. The original exploratory scripts were consolidated into these numbered scripts. Archived source scripts are retained in `scripts/KEEP_THESIS_DOCUMENTATION/` and `scripts/ARCHIVE_SUPERSEDED/` for traceability.

Before running, use `check_inputs.py` to verify the main paths and required columns. After running, use `summarize_outputs.py` to generate a compact table of key output counts.

Before running, also verify:
- Column names in your actual colon annotation input match what script 02 expects.
- The `rosetta_deployment_info.rds` file contains `pred_cols` and `factor_levels` that correspond to the output of script 02.
- The motif PFM table uses the same motif name format as the ROSETTA functional motif BED (e.g. `ZNF384_MA1125.2`).
- Chromosome names in the colon annotation input use integer codes (1–25), not `chr`-prefixed names — script 02 handles the conversion.
- The CRC variant BED has exactly the column structure expected by script 05 (see the `OVERLAP_COLUMNS` constant in the script for the expected 17-column format after bedtools intersect).
