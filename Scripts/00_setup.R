#!/usr/bin/env Rscript

# ============================================================
# 00_setup.R
#
# Shared setup for the CD marker isoform diversity analysis.
# This script defines packages, paths, global settings, and
# small helper functions used by downstream scripts.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(rtracklayer)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# PROJECT PATHS
# ============================================================
project_dir <- Sys.getenv("CD_MARKER_PROJECT_DIR", unset = getwd())

data_dir    <- file.path(project_dir, "data")
raw_dir     <- file.path(data_dir, "raw")
ref_dir     <- file.path(data_dir, "reference")
proc_dir    <- file.path(data_dir, "processed")
results_dir <- file.path(project_dir, "results")
fig_dir     <- file.path(project_dir, "figures")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# INPUT FILES
# ============================================================
# These are expected local filenames after downloading the
# external data/resources described in the README.

obj_path <- file.path(raw_dir, "combined_runs.rds")

gtf_path <- file.path(
  ref_dir,
  "Homo_sapiens.GRCh38.115.gtf.gz"
)

cd_marker_master_path <- file.path(
  ref_dir,
  "cd_marker_list.csv"
)

# ============================================================
# ANALYSIS SETTINGS
# ============================================================

ensembl_release <- 115L

min_cells_per_donor       <- 3L
min_donors                <- 2L
single_donor_rescue_cells <- 10L
donor_col                 <- "donor"

celltype_col <- "broad_celltype_final"

# ============================================================
# HELPER FUNCTIONS
# ============================================================

normalize_enst <- function(x) {
  stringr::str_extract(as.character(x), "ENST\\d+")
}

normalize_ensg <- function(x) {
  stringr::str_extract(as.character(x), "ENSG\\d+")
}

strip_version <- function(x) {
  sub("\\.\\d+$", "", x)
}

get_data_mat <- function(
    seu,
    assay = "Isoform",
    layer = "data",
    slot_fallback = "data"
) {
  if (!assay %in% names(seu@assays)) {
    return(NULL)
  }

  out <- tryCatch(
    Seurat::GetAssayData(seu, assay = assay, layer = layer),
    error = function(e) NULL
  )

  if (!is.null(out)) {
    return(out)
  }

  out <- tryCatch(
    Seurat::GetAssayData(seu, assay = assay, slot = slot_fallback),
    error = function(e) NULL
  )

  out
}

get_isoform_counts_matrix <- function(obj) {
  out <- tryCatch(
    SeuratObject::LayerData(obj, assay = "Isoform", layer = "counts"),
    error = function(e) NULL
  )

  if (!is.null(out)) {
    message("Using Isoform assay layer = counts")
    return(out)
  }

  out <- tryCatch(
    Seurat::GetAssayData(obj, assay = "Isoform", slot = "counts"),
    error = function(e) NULL
  )

  if (!is.null(out)) {
    message("Using Isoform assay slot = counts")
    return(out)
  }

  stop("Could not retrieve Isoform counts matrix.")
}

get_core_transcripts_by_donor <- function(
    obj,
    min_cells_per_donor = 3L,
    min_donors = 2L,
    single_donor_rescue_cells = 10L,
    donor_col = "donor"
) {
  iso_mat <- get_isoform_counts_matrix(obj)

  if (!donor_col %in% colnames(obj@meta.data)) {
    candidates <- c("sample", "patient", "orig.ident", "Donor", "Sample", "run")
    found <- candidates[candidates %in% colnames(obj@meta.data)]

    if (length(found) == 0) {
      stop(
        "No donor column found. Available metadata columns:\n",
        paste(colnames(obj@meta.data), collapse = ", ")
      )
    }

    donor_col <- found[1]
    message("Using '", donor_col, "' as donor column")
  }

  donors <- as.character(obj@meta.data[[donor_col]])
  unique_donors <- sort(unique(donors))

  message("Donors found: ", paste(unique_donors, collapse = ", "))
  message("Number of donors: ", length(unique_donors))

  donor_detection <- lapply(unique_donors, function(d) {
    cells_d <- which(donors == d)
    mat_d <- iso_mat[, cells_d, drop = FALSE]

    n_cells <- Matrix::rowSums(mat_d > 0)

    tibble::tibble(
      feature_id    = rownames(mat_d),
      transcript_id = normalize_enst(rownames(mat_d)),
      donor         = d,
      n_cells       = as.integer(n_cells)
    ) %>%
      dplyr::filter(!is.na(transcript_id), transcript_id != "")
  }) %>%
    dplyr::bind_rows()

  tx_donor_summary <- donor_detection %>%
    dplyr::mutate(detected = n_cells >= min_cells_per_donor) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(
      n_donors_detected   = sum(detected),
      donors_detected     = paste(sort(unique(donor[detected])), collapse = ","),
      max_cells_any_donor = max(n_cells),
      mean_cells_detected = ifelse(any(detected), mean(n_cells[detected]), NA_real_),
      .groups = "drop"
    )

  core_tx <- tx_donor_summary %>%
    dplyr::mutate(
      keep_core   = n_donors_detected >= min_donors,
      keep_rescue = max_cells_any_donor >= single_donor_rescue_cells,
      keep_final  = keep_core | keep_rescue,
      keep_reason = dplyr::case_when(
        keep_core & keep_rescue ~ "core_and_rescue",
        keep_core ~ "core_multi_donor",
        keep_rescue ~ "rescued_strong_single_donor",
        TRUE ~ "drop"
      )
    ) %>%
    dplyr::filter(keep_final)

  list(
    donor_detection  = donor_detection,
    tx_donor_summary = tx_donor_summary,
    core_tx          = core_tx
  )
}
