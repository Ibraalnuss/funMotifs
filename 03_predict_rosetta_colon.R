#!/usr/bin/env Rscript

# 03_predict_rosetta_colon.R
#
# Apply a trained ROSETTA rule model to prepared colon motif chunks.
#
# For each chunk TSV produced by 02_prepare_rosetta_colon_input.py, this script:
#   1. Loads the chunk and checks that all required predictor columns are present.
#   2. Coerces predictor columns to the factor levels used during training.
#   3. Runs predictClass() from R.ROSETTA to assign functional/non-functional labels.
#   4. Writes per-chunk full prediction files and functional-only files.
#   5. Appends functional motif rows to a combined BED file.
#   6. Writes a summary table with row counts and timing per chunk.
#
# Requires:
#   R.ROSETTA   installed in the R library
#   data.table  installed in the R library
#
# Example:
#   Rscript scripts/03_predict_rosetta_colon.R \
#     --chunk-dir data/interim/colon_rosetta_chunks \
#     --rules data/external/sig_rules_final_with_pretty_rules.rds \
#     --deployment data/external/rosetta_deployment_info.rds \
#     --out-dir results/tables/rosetta_chunk_predictions

suppressPackageStartupMessages({
  library(R.ROSETTA)
  library(data.table)
})

# ------------------------------------------------------------
# Small argument parser without extra dependencies
# ------------------------------------------------------------

get_arg <- function(args, key, default = NULL) {
  hit <- which(args == key)

  if (length(hit) == 0) {
    return(default)
  }

  if (hit == length(args)) {
    stop(paste("Missing value after", key))
  }

  args[hit + 1]
}

has_flag <- function(args, key) {
  key %in% args
}

args <- commandArgs(trailingOnly = TRUE)

if (has_flag(args, "--help") || length(args) == 0) {
  cat("
Usage:
  Rscript scripts/03_predict_rosetta_colon.R \\
    --chunk-dir <chunk_directory> \\
    --rules <sig_rules.rds> \\
    --deployment <rosetta_deployment_info.rds> \\
    --out-dir <output_directory>

Optional:
  --chunk-list <file_with_chunk_paths>
  --write-full TRUE/FALSE        default: TRUE
  --write-functional TRUE/FALSE  default: TRUE
  --write-bed TRUE/FALSE         default: TRUE
  --bed-name <filename>          default: colon_rosetta_functional_motifs.bed
  --summary-name <filename>      default: rosetta_prediction_summary.tsv

Outputs:
  Per chunk full prediction files, if --write-full TRUE
  Per chunk functional only files, if --write-functional TRUE
  One combined functional BED file, if --write-bed TRUE
  One summary table with row counts per chunk
")
  quit(status = 0)
}

chunk_dir <- get_arg(args, "--chunk-dir", NULL)
chunk_list <- get_arg(args, "--chunk-list", NULL)
rules_path <- get_arg(args, "--rules", "sig_rules_final_with_pretty_rules.rds")
deployment_path <- get_arg(args, "--deployment", "rosetta_deployment_info.rds")
out_dir <- get_arg(args, "--out-dir", "rosetta_predictions")

write_full <- tolower(get_arg(args, "--write-full", "TRUE")) == "true"
write_functional <- tolower(get_arg(args, "--write-functional", "TRUE")) == "true"
write_bed <- tolower(get_arg(args, "--write-bed", "TRUE")) == "true"

bed_name <- get_arg(args, "--bed-name", "colon_rosetta_functional_motifs.bed")
summary_name <- get_arg(args, "--summary-name", "rosetta_prediction_summary.tsv")

out_dir <- normalizePath(out_dir, mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rules_path)) {
  stop(paste("Rules file not found:", rules_path))
}

if (!file.exists(deployment_path)) {
  stop(paste("Deployment info file not found:", deployment_path))
}

# ------------------------------------------------------------
# Resolve chunk files
# ------------------------------------------------------------

if (!is.null(chunk_list)) {
  if (!file.exists(chunk_list)) {
    stop(paste("Chunk list file not found:", chunk_list))
  }

  chunk_files <- readLines(chunk_list)
  chunk_files <- chunk_files[nchar(chunk_files) > 0]

} else {
  if (is.null(chunk_dir)) {
    stop("Provide either --chunk-dir or --chunk-list.")
  }

  if (!dir.exists(chunk_dir)) {
    stop(paste("Chunk directory not found:", chunk_dir))
  }

  chunk_files <- list.files(
    chunk_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  )
}

chunk_files <- sort(chunk_files)

if (length(chunk_files) == 0) {
  stop("No chunk files found.")
}

cat("Chunks found:", length(chunk_files), "\n")

# ------------------------------------------------------------
# Load ROSETTA objects
# ------------------------------------------------------------

cat("Loading rules:", rules_path, "\n")
sig_rules <- readRDS(rules_path)

cat("Loading deployment info:", deployment_path, "\n")
deployment_info <- readRDS(deployment_path)

pred_cols <- deployment_info$pred_cols
factor_levels <- deployment_info$factor_levels

if (is.null(pred_cols) || is.null(factor_levels)) {
  stop("Deployment info must contain pred_cols and factor_levels.")
}

meta_cols <- c("row_id", "mid", "chr", "motifstart", "motifend", "name", "strand")

# ------------------------------------------------------------
# Prepare combined output files
# ------------------------------------------------------------

combined_bed_path <- file.path(out_dir, bed_name)
summary_path <- file.path(out_dir, summary_name)

if (file.exists(combined_bed_path)) {
  file.remove(combined_bed_path)
}

summary_rows <- list()

# ------------------------------------------------------------
# Prediction helper
# ------------------------------------------------------------

predict_one_chunk <- function(chunk_path, chunk_index) {
  cat("\nProcessing chunk", chunk_index, ":", chunk_path, "\n")

  t0 <- Sys.time()

  df <- fread(
    chunk_path,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    showProgress = FALSE
  )

  t_read <- Sys.time()

  missing_meta <- setdiff(meta_cols, colnames(df))
  if (length(missing_meta) > 0) {
    stop(paste(
      "Missing metadata columns in",
      chunk_path,
      ":",
      paste(missing_meta, collapse = ", ")
    ))
  }

  missing_pred <- setdiff(pred_cols, colnames(df))
  if (length(missing_pred) > 0) {
    stop(paste(
      "Missing predictor columns in",
      chunk_path,
      ":",
      paste(missing_pred, collapse = ", ")
    ))
  }

  x <- df[, pred_cols, drop = FALSE]

  for (nm in pred_cols) {
    x[[nm]] <- factor(as.character(x[[nm]]), levels = factor_levels[[nm]])
  }

  invalid_counts <- sapply(x, function(col) sum(is.na(col)))

  if (any(invalid_counts > 0)) {
    print(invalid_counts[invalid_counts > 0])
    stop(paste("Factor coercion produced NA values in", chunk_path))
  }

  t_factor <- Sys.time()

  pred_obj <- predictClass(
    dt = x,
    rules = sig_rules,
    discrete = TRUE,
    validate = FALSE
  )

  t_pred <- Sys.time()

  pred_df <- pred_obj$out

  if (nrow(pred_df) != nrow(df)) {
    stop(paste("Prediction row count mismatch in", chunk_path))
  }

  result <- cbind(df[, meta_cols, drop = FALSE], pred_df)
  functional <- result[result$predictedClass == "1", , drop = FALSE]

  base_name <- tools::file_path_sans_ext(basename(chunk_path))

  full_out <- file.path(out_dir, paste0(base_name, "_predictions.tsv"))
  functional_out <- file.path(out_dir, paste0(base_name, "_functional_only.tsv"))

  if (write_full) {
    fwrite(
      result,
      file = full_out,
      sep = "\t",
      quote = FALSE
    )
  }

  if (write_functional) {
    fwrite(
      functional,
      file = functional_out,
      sep = "\t",
      quote = FALSE
    )
  }

  if (write_bed && nrow(functional) > 0) {
    bed <- functional[, c("chr", "motifstart", "motifend", "name", "strand", "mid")]

    fwrite(
      bed,
      file = combined_bed_path,
      sep = "\t",
      quote = FALSE,
      col.names = FALSE,
      append = file.exists(combined_bed_path)
    )
  }

  t_end <- Sys.time()

  predicted_counts <- table(result$predictedClass)

  cat("Rows:", format(nrow(result), big.mark = ","), "\n")
  cat("Functional:", format(nrow(functional), big.mark = ","), "\n")
  cat("Read sec:", round(as.numeric(difftime(t_read, t0, units = "secs")), 2), "\n")
  cat("Factor sec:", round(as.numeric(difftime(t_factor, t_read, units = "secs")), 2), "\n")
  cat("Predict sec:", round(as.numeric(difftime(t_pred, t_factor, units = "secs")), 2), "\n")
  cat("Total sec:", round(as.numeric(difftime(t_end, t0, units = "secs")), 2), "\n")

  summary <- data.frame(
    chunk_index = chunk_index,
    chunk_file = chunk_path,
    rows = nrow(result),
    functional_rows = nrow(functional),
    class_0 = ifelse("0" %in% names(predicted_counts), as.integer(predicted_counts[["0"]]), 0),
    class_1 = ifelse("1" %in% names(predicted_counts), as.integer(predicted_counts[["1"]]), 0),
    read_sec = as.numeric(difftime(t_read, t0, units = "secs")),
    factor_sec = as.numeric(difftime(t_factor, t_read, units = "secs")),
    predict_sec = as.numeric(difftime(t_pred, t_factor, units = "secs")),
    total_sec = as.numeric(difftime(t_end, t0, units = "secs"))
  )

  rm(df, x, pred_obj, pred_df, result, functional)
  gc()

  summary
}

# ------------------------------------------------------------
# Run all chunks
# ------------------------------------------------------------

for (i in seq_along(chunk_files)) {
  summary_rows[[i]] <- predict_one_chunk(chunk_files[[i]], i)
}

summary_df <- rbindlist(summary_rows, fill = TRUE)

fwrite(
  summary_df,
  file = summary_path,
  sep = "\t",
  quote = FALSE
)

cat("\nDone.\n")
cat("Chunks processed:", nrow(summary_df), "\n")
cat("Total rows scored:", format(sum(summary_df$rows), big.mark = ","), "\n")
cat("Total functional motifs:", format(sum(summary_df$functional_rows), big.mark = ","), "\n")
cat("Summary:", summary_path, "\n")

if (write_bed) {
  cat("Combined functional BED:", combined_bed_path, "\n")
}
