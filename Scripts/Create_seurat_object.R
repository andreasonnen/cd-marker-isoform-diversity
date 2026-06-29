#!/usr/bin/env Rscript

# ============================================================
# Build the Seurat object from Ensembl115 kb/kallisto quant-tcc outputs.
#
# Input per sample:
#   - matrix.abundance.mtx
#   - matrix.abundance.tpm.mtx
#   - matrix.abundance.gene.mtx
#   - matrix.abundance.gene.tpm.mtx
#   - transcripts.txt
#   - genes.txt
#   - barcode file from kb/counts_unfiltered
#
# Output:
#   - Seurat object with:
#       Isoform assay: transcript-level counts + TPM
#       Gene assay: gene-level counts + TPM + normalized data
#       gene annotation from Ensembl 115 GTF
#       cell-level QC metrics
#       broad cell-type labels
#       PCA, Harmony, UMAP and clusters
#
# Example sample sheet row:
#   r6d6,/path/to/Isoform_counts_filtered,/path/to/cells_x_tcc.barcodes.txt,run6,donor6,r6d6,TRUE
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(harmony)
  library(patchwork)
  library(purrr)
})

options(Seurat.object.assay.version = "v5")
set.seed(1)

# ============================================================
# CONFIG
# ============================================================

# These can be edited directly or supplied as environment variables.
sample_sheet <- Sys.getenv(
  "CD_SAMPLE_SHEET",
  unset = "../Data/kb_quant_sample_sheet.csv"
)

barcode_annotation_file <- Sys.getenv(
  "CD_BARCODE_ANNOTATION",
  unset = "../Data/PBMCs.allruns.barcode_annotation.txt"
)

gtf_path <- Sys.getenv(
  "CD_ENSEMBL115_GTF",
  unset = "../Data/Homo_sapiens.GRCh38.115.clean.gtf"
)

out_dir <- Sys.getenv(
  "CD_OUT_DIR",
  unset = "../Results/preprocessing"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(out_dir, "combined_runs_preprocessed_ensembl115.rds")
out_metadata <- file.path(out_dir, "combined_runs_cell_metadata.csv")
out_hvg <- file.path(out_dir, "combined_runs_HVGs_after_TRIG_removal.txt")

isoform_assay <- "Isoform"
gene_assay <- "Gene"

# Cell type column in the deposited Smart-seq3xpress barcode annotation.
celltype_col <- "celltype_lvl2_inex_10khvg_reads_res08_new"

# Final QC thresholds.
qc_min_features <- 500L
qc_max_features <- 7500L
qc_max_gene_counts <- 150000
qc_max_percent_mt <- 8

# Gene-level dimensionality reduction settings.
n_pcs <- 20
cluster_resolution <- 0.2
harmony_vars <- c("donor")

# ============================================================
# SAMPLE SHEET TEMPLATE
# ============================================================

if (!file.exists(sample_sheet)) {
  template <- tibble::tibble(
    sample_id = c("r6d6", "r6d7", "r6d8"),
    isoform_dir = c(
      "/path/to/donor6/kb_ss3xpress_ensembl115_filtered/Isoform_counts_filtered",
      "/path/to/donor7/kb_ss3xpress_ensembl115_filtered/Isoform_counts_filtered",
      "/path/to/donor8/kb_ss3xpress_ensembl115_filtered/Isoform_counts_filtered"
    ),
    barcode_path = c(
      "/path/to/donor6/kb_ss3xpress_ensembl115_filtered/counts_unfiltered/cells_x_tcc.barcodes.txt",
      "/path/to/donor7/kb_ss3xpress_ensembl115_filtered/counts_unfiltered/cells_x_tcc.barcodes.txt",
      "/path/to/donor8/kb_ss3xpress_ensembl115_filtered/counts_unfiltered/cells_x_tcc.barcodes.txt"
    ),
    run = c("run6", "run6", "run6"),
    donor = c("donor6", "donor7", "donor8"),
    id_prefix = c("r6d6", "r6d7", "r6d8"),
    include = c(TRUE, TRUE, TRUE)
  )

  dir.create(dirname(sample_sheet), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(template, sample_sheet)

  stop(
    "Sample sheet not found. A template was written to:\n",
    sample_sheet,
    "\nEdit this file with the correct paths and rerun the script."
  )
}

# ============================================================
# HELPERS
# ============================================================

read_kb_pair <- function(count_path, tpm_path, feature_path, barcode_path) {
  counts_raw <- Matrix::readMM(count_path)
  tpm_raw <- Matrix::readMM(tpm_path)

  features <- make.unique(readLines(feature_path))
  barcodes <- readLines(barcode_path)

  if (!all(dim(counts_raw) == dim(tpm_raw))) {
    stop("Counts and TPM matrices do not have the same dimensions.")
  }

  # Case A: rows = cells, columns = features
  if (nrow(counts_raw) == length(barcodes) && ncol(counts_raw) == length(features)) {
    counts <- Matrix::t(counts_raw)
    tpm <- Matrix::t(tpm_raw)

    # Case B: rows = features, columns = cells
  } else if (nrow(counts_raw) == length(features) && ncol(counts_raw) == length(barcodes)) {
    counts <- counts_raw
    tpm <- tpm_raw

  } else {
    stop(
      paste0(
        "Matrix dimensions do not match feature/barcode lengths.\n",
        "counts dim: ", paste(dim(counts_raw), collapse = " x "), "\n",
        "TPM dim: ", paste(dim(tpm_raw), collapse = " x "), "\n",
        "n features: ", length(features), "\n",
        "n barcodes: ", length(barcodes), "\n",
        "Check that the barcode file belongs to this quantification."
      )
    )
  }

  rownames(counts) <- features
  colnames(counts) <- barcodes

  rownames(tpm) <- features
  colnames(tpm) <- barcodes

  counts <- methods::as(counts, "CsparseMatrix")
  tpm <- methods::as(tpm, "CsparseMatrix")

  list(counts = counts, tpm = tpm)
}

load_run_with_gene_and_isoform <- function(
    isoform_dir,
    barcode_path,
    run_name,
    donor_name,
    id_prefix,
    isoform_assay = "Isoform",
    gene_assay = "Gene"
) {
  message("Loading sample: ", id_prefix)

  iso <- read_kb_pair(
    count_path = file.path(isoform_dir, "matrix.abundance.mtx"),
    tpm_path = file.path(isoform_dir, "matrix.abundance.tpm.mtx"),
    feature_path = file.path(isoform_dir, "transcripts.txt"),
    barcode_path = barcode_path
  )

  gene <- read_kb_pair(
    count_path = file.path(isoform_dir, "matrix.abundance.gene.mtx"),
    tpm_path = file.path(isoform_dir, "matrix.abundance.gene.tpm.mtx"),
    feature_path = file.path(isoform_dir, "genes.txt"),
    barcode_path = barcode_path
  )

  if (!identical(colnames(iso$counts), colnames(gene$counts))) {
    stop("Isoform and gene matrices do not have identical cell columns/order.")
  }

  obj <- CreateSeuratObject(
    counts = iso$counts,
    assay = isoform_assay,
    project = "SS3_PBMC",
    min.cells = 0,
    min.features = 0
  )

  # Preserve original barcode before renaming cells.
  obj$barcode <- colnames(obj)

  # Add transcript TPM as an explicit layer.
  obj[[isoform_assay]]$tpm <- iso$tpm

  # Add gene assay with count and TPM layers.
  gene_assay_obj <- CreateAssay5Object(
    counts = gene$counts,
    min.cells = 0,
    min.features = 0
  )

  obj[[gene_assay]] <- gene_assay_obj
  obj[[gene_assay]]$tpm <- gene$tpm

  obj$run <- run_name
  obj$donor <- donor_name
  obj$sample_id <- paste(run_name, donor_name, sep = "_")

  new_cell_names <- paste0(id_prefix, "_", colnames(obj))
  obj <- RenameCells(obj, new.names = new_cell_names)

  DefaultAssay(obj) <- gene_assay
  obj
}

join_layer_if_possible <- function(obj, assay, layer) {
  DefaultAssay(obj) <- assay

  obj[[assay]] <- tryCatch(
    {
      JoinLayers(obj[[assay]], layers = layer)
    },
    error = function(e) {
      message("Skipping JoinLayers for ", assay, " / ", layer, ": ", e$message)
      obj[[assay]]
    }
  )

  obj
}

add_assay_count_metrics <- function(obj, assay) {
  counts <- LayerData(obj[[assay]], layer = "counts")

  obj[[paste0("nCount_", assay)]] <- Matrix::colSums(counts)
  obj[[paste0("nFeature_", assay)]] <- Matrix::colSums(counts > 0)

  obj
}

extract_gtf_attr <- function(x, key) {
  pat <- paste0(key, ' "([^"]+)"')
  ifelse(
    grepl(pat, x),
    sub(paste0(".*", pat, ".*"), "\\1", x),
    NA_character_
  )
}

read_gtf_gene_map <- function(gtf_path) {
  cmd <- if (grepl("\\.gz$", gtf_path)) {
    paste("zcat", shQuote(gtf_path), "| grep -v '^#'")
  } else {
    paste("grep -v '^#'", shQuote(gtf_path))
  }

  gtf <- data.table::fread(
    cmd = cmd,
    sep = "\t",
    header = FALSE,
    quote = "",
    data.table = FALSE
  )

  colnames(gtf) <- c(
    "seqname", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute"
  )

  genes <- gtf %>%
    dplyr::filter(feature == "gene") %>%
    dplyr::transmute(
      gene_id = sub("\\..*$", "", extract_gtf_attr(attribute, "gene_id")),
      gene_name = extract_gtf_attr(attribute, "gene_name"),
      gene_biotype = ifelse(
        grepl('gene_biotype "', attribute),
        extract_gtf_attr(attribute, "gene_biotype"),
        extract_gtf_attr(attribute, "gene_type")
      )
    ) %>%
    dplyr::distinct()

  genes
}

annotate_gene_assay <- function(obj, gene_annot, gene_assay = "Gene") {
  DefaultAssay(obj) <- gene_assay

  feat_md <- obj[[gene_assay]][[]]
  current_rows <- rownames(obj[[gene_assay]])

  if ("ensembl_gene_id" %in% colnames(feat_md)) {
    gene_ids <- feat_md$ensembl_gene_id
    bad <- is.na(gene_ids) | gene_ids == ""
    gene_ids[bad] <- current_rows[bad]
  } else {
    gene_ids <- current_rows
  }

  gene_ids <- sub("\\..*$", "", gene_ids)

  bad2 <- is.na(gene_ids) | gene_ids == ""
  gene_ids[bad2] <- current_rows[bad2]
  gene_ids <- sub("\\..*$", "", gene_ids)

  gene_name_vec <- gene_annot$gene_name[match(gene_ids, gene_annot$gene_id)]
  gene_biotype_vec <- gene_annot$gene_biotype[match(gene_ids, gene_annot$gene_id)]

  name_counts <- table(gene_name_vec[!is.na(gene_name_vec) & gene_name_vec != ""])

  dup_flag <- ifelse(
    !is.na(gene_name_vec) & gene_name_vec != "",
    name_counts[gene_name_vec] > 1,
    FALSE
  )

  new_gene_names <- ifelse(
    is.na(gene_name_vec) | gene_name_vec == "",
    gene_ids,
    ifelse(
      dup_flag,
      paste0(gene_name_vec, "|", gene_ids),
      gene_name_vec
    )
  )

  bad3 <- is.na(new_gene_names) | new_gene_names == ""
  new_gene_names[bad3] <- gene_ids[bad3]

  bad4 <- is.na(new_gene_names) | new_gene_names == ""
  new_gene_names[bad4] <- current_rows[bad4]

  new_gene_names <- make.unique(new_gene_names)

  dimnames(obj[[gene_assay]]) <- list(new_gene_names, colnames(obj))

  gene_lookup <- data.frame(
    gene_name_display = new_gene_names,
    ensembl_gene_id = gene_ids,
    gene_name_raw = gene_name_vec,
    gene_biotype = gene_biotype_vec,
    stringsAsFactors = FALSE,
    row.names = new_gene_names
  )

  obj[[gene_assay]] <- AddMetaData(
    object = obj[[gene_assay]],
    metadata = gene_lookup
  )

  obj
}

add_gene_qc_metrics <- function(obj, gene_assay = "Gene") {
  DefaultAssay(obj) <- gene_assay

  feat_md <- obj[[gene_assay]][[]]

  mt_features <- rownames(obj[[gene_assay]])[
    !is.na(feat_md$gene_name_raw) &
      grepl("^MT-", feat_md$gene_name_raw)
  ]

  ribo_features <- rownames(obj[[gene_assay]])[
    !is.na(feat_md$gene_name_raw) &
      grepl("^RPS([0-9]+[A-Z]?|A)$|^RPL([0-9]+[A-Z]?|P[0-2])$", feat_md$gene_name_raw) &
      feat_md$gene_biotype == "protein_coding"
  ]

  message("Mitochondrial features: ", length(mt_features))
  message("Ribosomal features: ", length(ribo_features))

  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    features = mt_features,
    assay = gene_assay
  )

  obj[["percent.ribo"]] <- PercentageFeatureSet(
    obj,
    features = ribo_features,
    assay = gene_assay
  )

  obj
}

add_barcode_metadata <- function(obj, anno_df, celltype_col) {
  if (!"barcode" %in% colnames(obj@meta.data)) {
    stop("Object metadata does not contain original barcode column.")
  }

  obj$barcode_run_donor <- paste(
    obj$barcode,
    obj$run,
    obj$donor,
    sep = "_"
  )

  needed_cols <- c("barcode_run_donor", "QC_status", celltype_col)
  missing_cols <- setdiff(needed_cols, colnames(anno_df))

  if (length(missing_cols) > 0) {
    stop("Missing columns in barcode annotation: ", paste(missing_cols, collapse = ", "))
  }

  lookup <- anno_df[, needed_cols, drop = FALSE]
  rownames(lookup) <- lookup$barcode_run_donor

  obj <- AddMetaData(
    obj,
    metadata = lookup[
      obj$barcode_run_donor,
      c("QC_status", celltype_col),
      drop = FALSE
    ]
  )

  obj$celltype <- obj[[celltype_col]][, 1]
  obj[[celltype_col]] <- NULL

  obj
}

collapse_to_broad_celltypes <- function(obj) {
  obj$broad_celltype_final <- dplyr::case_when(
    obj$celltype %in% c("Naive B", "Memory B") ~ "B_cell",

    obj$celltype %in% c(
      "Naive CD4+ T",
      "CD4+ TCM/TEM",
      "Naive/CM CD4+ T_1",
      "Naive/CM CD4+ T_2",
      "Clonal CD4+ T",
      "Naive/CM CD4/CD8+ T"
    ) ~ "CD4_T",

    obj$celltype %in% c(
      "Naive CD8+ T",
      "CD8+ TCM/TEM",
      "Clonal CD8+ T",
      "CD8+ TEM/CD4+ CTL"
    ) ~ "CD8_T",

    obj$celltype == "Tregs" ~ "Treg",

    obj$celltype %in% c("NK_1", "NK_2") ~ "NK",
    obj$celltype == "NK/ILC" ~ "NK/ILC",

    obj$celltype == "MAIT" ~ "MAIT",
    obj$celltype %in% c("gdT_1", "gdT_2") ~ "gdT",

    obj$celltype %in% c("CD14+ Mono", "CD16+ Mono") ~ "Monocyte",
    obj$celltype %in% c("cDC_1", "cDC_2") ~ "DC",
    obj$celltype == "pDC" ~ "pDC",

    obj$celltype == "HSPCs" ~ "HSPC",
    obj$celltype == "Platelets" ~ "Platelet",

    obj$celltype == "Proliferating cells" ~ "Proliferating",

    TRUE ~ NA_character_
  )

  obj$broad_celltype_final <- factor(
    obj$broad_celltype_final,
    levels = c(
      "CD4_T", "CD8_T", "Treg", "B_cell", "NK",
      "DC", "Monocyte", "pDC", "NK/ILC", "gdT",
      "MAIT", "HSPC", "Platelet", "Proliferating"
    )
  )

  obj
}

save_qc_plots <- function(obj, out_dir, gene_assay = "Gene") {
  DefaultAssay(obj) <- gene_assay

  p_vln <- VlnPlot(
    obj,
    features = c("nFeature_Gene", "nCount_Gene", "percent.mt", "percent.ribo"),
    group.by = "sample_id",
    ncol = 4,
    pt.size = 0
  )

  ggsave(
    file.path(out_dir, "QC_violin_gene_metrics_by_sample.pdf"),
    p_vln,
    width = 14,
    height = 5,
    device = cairo_pdf,
    bg = "white"
  )

  p_scatter1 <- FeatureScatter(
    obj,
    feature1 = "nCount_Gene",
    feature2 = "nFeature_Gene"
  )

  p_scatter2 <- FeatureScatter(
    obj,
    feature1 = "nCount_Gene",
    feature2 = "percent.mt"
  )

  p_scatter3 <- FeatureScatter(
    obj,
    feature1 = "nFeature_Gene",
    feature2 = "percent.mt"
  )

  p_scatter <- p_scatter1 | p_scatter2 | p_scatter3

  ggsave(
    file.path(out_dir, "QC_scatter_gene_metrics.pdf"),
    p_scatter,
    width = 14,
    height = 4.5,
    device = cairo_pdf,
    bg = "white"
  )
}

# ============================================================
# 1) LOAD SAMPLE SHEET
# ============================================================

samples <- readr::read_csv(sample_sheet, show_col_types = FALSE) %>%
  dplyr::mutate(include = as.logical(include)) %>%
  dplyr::filter(include)

required_sample_cols <- c(
  "sample_id", "isoform_dir", "barcode_path",
  "run", "donor", "id_prefix", "include"
)

missing_sample_cols <- setdiff(required_sample_cols, colnames(samples))

if (length(missing_sample_cols) > 0) {
  stop("Missing columns in sample sheet: ", paste(missing_sample_cols, collapse = ", "))
}

message("Samples to load:")
print(samples %>% dplyr::select(sample_id, run, donor, id_prefix), n = Inf)

# ============================================================
# 2) LOAD PER-SAMPLE OBJECTS
# ============================================================

obj_list <- purrr::pmap(
  list(
    isoform_dir = samples$isoform_dir,
    barcode_path = samples$barcode_path,
    run_name = samples$run,
    donor_name = samples$donor,
    id_prefix = samples$id_prefix
  ),
  function(isoform_dir, barcode_path, run_name, donor_name, id_prefix) {
    load_run_with_gene_and_isoform(
      isoform_dir = isoform_dir,
      barcode_path = barcode_path,
      run_name = run_name,
      donor_name = donor_name,
      id_prefix = id_prefix,
      isoform_assay = isoform_assay,
      gene_assay = gene_assay
    )
  }
)

names(obj_list) <- samples$sample_id

# ============================================================
# 3) MERGE OBJECTS BEFORE GENE-SYMBOL RENAMING
# ============================================================

message("Merging samples...")

combined <- Reduce(
  f = function(x, y) {
    merge(x = x, y = y, project = "SS3_PBMC")
  },
  x = obj_list
)

# ============================================================
# 4) IMPORT CELL-LEVEL METADATA AND KEEP QCPASS CELLS
# ============================================================

message("Loading barcode annotation metadata...")

anno_all <- read.delim(
  barcode_annotation_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

anno_all$barcode_run_donor <- paste(
  anno_all$barcode,
  anno_all$run,
  anno_all$donor,
  sep = "_"
)

combined <- add_barcode_metadata(
  obj = combined,
  anno_df = anno_all,
  celltype_col = celltype_col
)

message("Cells before QCpass filter: ", ncol(combined))

combined <- subset(
  combined,
  subset = !is.na(QC_status) & QC_status == "QCpass"
)

message("Cells after QCpass filter: ", ncol(combined))

# ============================================================
# 5) JOIN SEURAT V5 LAYERS AFTER MERGE/SUBSETTING
# ============================================================

combined <- join_layer_if_possible(combined, assay = isoform_assay, layer = "counts")
combined <- join_layer_if_possible(combined, assay = isoform_assay, layer = "tpm")

combined <- join_layer_if_possible(combined, assay = gene_assay, layer = "counts")
combined <- join_layer_if_possible(combined, assay = gene_assay, layer = "tpm")

combined <- add_assay_count_metrics(combined, assay = isoform_assay)
combined <- add_assay_count_metrics(combined, assay = gene_assay)

# ============================================================
# 6) ANNOTATE GENE ASSAY FROM ENSEMBL 115 GTF
# ============================================================

message("Reading GTF gene annotation...")
gene_annot <- read_gtf_gene_map(gtf_path)

combined <- annotate_gene_assay(
  obj = combined,
  gene_annot = gene_annot,
  gene_assay = gene_assay
)

# Sanity checks after gene annotation.
feat_md <- combined[[gene_assay]][[]]
gene_rows <- rownames(combined[[gene_assay]])

message("Gene assay features: ", length(gene_rows))
message("Feature metadata rows: ", nrow(feat_md))
message("Gene rownames identical to feature metadata rownames: ", identical(gene_rows, rownames(feat_md)))
message("Duplicated gene rownames: ", sum(duplicated(gene_rows)))
message("Rows still starting with ENSG: ", sum(grepl("^ENSG", gene_rows)))
message("Rows starting with MT-: ", sum(grepl("^MT-", gene_rows)))

# ============================================================
# 7) ADD QC METRICS AND FILTER CELLS
# ============================================================

combined <- add_gene_qc_metrics(
  obj = combined,
  gene_assay = gene_assay
)

message("QC summaries before additional filtering:")
print(summary(combined$nFeature_Gene))
print(summary(combined$nCount_Gene))
print(summary(combined$percent.mt))
print(summary(combined$percent.ribo))

save_qc_plots(
  obj = combined,
  out_dir = out_dir,
  gene_assay = gene_assay
)

message("Cells before additional QC filter: ", ncol(combined))

combined <- subset(
  combined,
  subset =
    nFeature_Gene >= qc_min_features &
    nFeature_Gene <= qc_max_features &
    nCount_Gene <= qc_max_gene_counts &
    percent.mt < qc_max_percent_mt
)

message("Cells after additional QC filter: ", ncol(combined))

message("QC summaries after filtering:")
print(summary(combined$nFeature_Gene))
print(summary(combined$nCount_Gene))
print(summary(combined$percent.mt))
print(summary(combined$percent.ribo))

# Recompute assay count metrics after filtering.
combined <- add_assay_count_metrics(combined, assay = isoform_assay)
combined <- add_assay_count_metrics(combined, assay = gene_assay)

combined$depth_gene_counts <- combined$nCount_Gene
combined$detected_genes <- combined$nFeature_Gene

combined$depth_isoform_counts <- combined$nCount_Isoform
combined$detected_isoforms <- combined$nFeature_Isoform

# ============================================================
# 8) GENE-LEVEL NORMALIZATION, HVG SELECTION, PCA
# ============================================================

DefaultAssay(combined) <- gene_assay

combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, verbose = FALSE)

# Remove TCR and immunoglobulin genes from the HVG set so that
# clonotype-like signal does not dominate the reduced space.
tr_ig_genes <- grep(
  pattern = paste(
    c(
      "^TRAV", "^TRAJ", "^TRBV", "^TRBJ",
      "^TRDV", "^TRDJ", "^TRGV", "^TRGJ",
      "^IGHV", "^IGHD", "^IGHJ",
      "^IGKV", "^IGKJ",
      "^IGLV", "^IGLJ"
    ),
    collapse = "|"
  ),
  x = VariableFeatures(combined),
  value = TRUE
)

VariableFeatures(combined) <- setdiff(
  VariableFeatures(combined),
  tr_ig_genes
)

message("Variable features after TR/IG removal: ", length(VariableFeatures(combined)))

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, verbose = FALSE)

stdev <- Stdev(combined, reduction = "pca")
var_explained <- (stdev^2) / sum(stdev^2) * 100

message("Cumulative variance explained by first ", n_pcs, " PCs: ",
        round(sum(var_explained[seq_len(n_pcs)]), 2), "%")

# ============================================================
# 9) HARMONY, NEIGHBOURS, CLUSTERS, UMAP
# ============================================================

pca_mat <- Embeddings(combined, "pca")[, seq_len(n_pcs), drop = FALSE]
md <- combined@meta.data[rownames(pca_mat), , drop = FALSE]

message("Running Harmony with variables: ", paste(harmony_vars, collapse = ", "))

harmony_mat <- harmony::HarmonyMatrix(
  data_mat = pca_mat,
  meta_data = md,
  vars_use = harmony_vars,
  do_pca = FALSE,
  plot_convergence = FALSE
)

combined[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_mat,
  key = "harmony_",
  assay = DefaultAssay(combined)
)

combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = seq_len(n_pcs),
  verbose = FALSE
)

combined <- FindClusters(
  combined,
  resolution = cluster_resolution,
  verbose = FALSE
)

combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = seq_len(n_pcs),
  verbose = FALSE
)

# ============================================================
# 10) BROAD CELL-TYPE LABELS
# ============================================================

combined <- collapse_to_broad_celltypes(combined)

message("Broad cell type counts:")
print(table(combined$broad_celltype_final, useNA = "ifany"))

# Save simple annotation plots.
p_celltype <- DimPlot(
  combined,
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
)

p_broad <- DimPlot(
  combined,
  group.by = "broad_celltype_final",
  label = TRUE,
  repel = TRUE
)

p_sample <- DimPlot(
  combined,
  group.by = "sample_id"
)

ggsave(
  file.path(out_dir, "UMAP_celltype_original.pdf"),
  p_celltype,
  width = 8,
  height = 6,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(out_dir, "UMAP_broad_celltype_final.pdf"),
  p_broad,
  width = 8,
  height = 6,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(out_dir, "UMAP_sample_id.pdf"),
  p_sample,
  width = 8,
  height = 6,
  device = cairo_pdf,
  bg = "white"
)

# Marker check plot for broad labels.
DefaultAssay(combined) <- gene_assay

markers_to_check <- c(
  # T cells
  "CD3D", "CD3E", "TRAC",
  "IL7R", "LTB", "CCR7", "TCF7", "LEF1",
  "CD8A", "CD8B", "CCL5", "PRF1", "GZMB", "NKG7",
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",

  # B cells
  "MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA",

  # Monocytes/DC/pDC
  "LYZ", "FCN1", "S100A8", "S100A9", "LST1",
  "FCGR3A", "MS4A7",
  "FCER1A", "CST3", "CLEC10A", "CD1C",
  "CLEC4C", "IL3RA", "TCF4", "JCHAIN",

  # NK/gdT/MAIT
  "GNLY", "KLRD1", "TYROBP", "XCL1", "XCL2",
  "TRDC", "TRGC1", "TRGC2",
  "KLRB1", "SLC4A10", "TRAV1-2",

  # HSPC/platelet/proliferating
  "AVP", "GATA2", "SOX4",
  "PPBP", "PF4",
  "MKI67", "TOP2A", "STMN1"
)

markers_to_check <- unique(markers_to_check)
markers_to_check <- markers_to_check[markers_to_check %in% rownames(combined[[gene_assay]])]

p_dot <- DotPlot(
  combined,
  features = markers_to_check,
  group.by = "broad_celltype_final",
  assay = gene_assay
) +
  RotatedAxis()

ggsave(
  file.path(out_dir, "DotPlot_broad_celltype_marker_check.pdf"),
  p_dot,
  width = 16,
  height = 6.5,
  device = cairo_pdf,
  bg = "white"
)

# ============================================================
# 11) SAVE OUTPUTS
# ============================================================

message("Saving object: ", out_rds)
saveRDS(combined, out_rds)

readr::write_csv(
  combined@meta.data %>%
    tibble::rownames_to_column("cell_id"),
  out_metadata
)

readr::write_lines(
  VariableFeatures(combined),
  out_hvg
)

message("Done.")
message("Final object:")
print(combined)
