#!/usr/bin/env Rscript

# Final Figures manuscript. CD markers paper

# ============================================================
# FIGURE 1: CD marker transcript complexity, coding potential,
# within-population coding transcript heterogeneity, and
# pairwise isoform switching
#
# Panels:
# A = transcript complexity per CD marker gene
# B = cell-type-specific coding potential
# C = strict within-population JS boxplot
# D = pairwise isoform switching triangle
#
# Saves:
# - Combined figure as PDF + PNG
# - Each panel separately as PDF + PNG
# - Panel B diagnostic tables
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(tibble)
  library(scales)
})

# ============================================================
# PATHS
# ============================================================

fig_dir <- "/home/projects/dtu_00062/people/andson/SS3_PBMC/CD_markers_paper/Ensembl115/Figures_Manuscript"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fig1_universe_file <- file.path(
  fig_dir,
  "Figure1_CDmarker_transcript_universe_hybrid_min5cells_min2donors_rescue10_fullGTF.tsv"
)

seurat_file_1 <- "/home/projects/dtu_00062/people/andson/SS3_PBMC/CD_markers_paper/Ensembl115/Data/combined_runs.rds"
seurat_file_2 <- "/home/projects/dtu_00062/people/andson/SS3_PBMC/Data/combined_runs.rds"
seurat_file <- if (file.exists(seurat_file_1)) seurat_file_1 else seurat_file_2

cd_marker_master_file <- "/home/projects/dtu_00062/people/andson/SS3_PBMC/CD_markers_paper/Ensembl115/Data/cd_marker_list_clean_ensembl115_updated.csv"

js_file <- file.path(
  fig_dir,
  "Figure1C_within_population_coding_transcript_JS_distances.csv"
)

diag_file <- file.path(
  fig_dir,
  "Figure1C_JS_sparsity_diagnostic.csv"
)

isar_base_dir <- "/home/projects/dtu_00062/people/andson/SS3_PBMC/CD_markers_paper/Ensembl115/Results/Transcript_switching_analysis/ISAR_results_two_pass_broad_celltype_final"
pass_pc_dir <- file.path(isar_base_dir, "protein_coding_only")

gene_calls_pc_file <- file.path(
  pass_pc_dir,
  "protein_coding_only_all_top_switch_genes.rds"
)

usable_ct_pc_file <- file.path(
  pass_pc_dir,
  "protein_coding_only_usable_celltypes.rds"
)

# ============================================================
# OUTPUT FILES
# ============================================================

outfile_combined_pdf <- file.path(fig_dir, "Figure1_ABCD_final.pdf")
outfile_combined_png <- file.path(fig_dir, "Figure1_ABCD_final.png")

# ============================================================
# SETTINGS
# ============================================================

isoform_assay <- "Isoform"
signal_layer  <- "counts"

celltype_var <- "broad_celltype_final"

celltype_map <- c(
  "CD4_T"    = "CD4 T",
  "CD8_T"    = "CD8 T",
  "B_cell"   = "B cell",
  "NK"       = "NK",
  "Treg"     = "Treg",
  "pDC"      = "pDC",
  "DC"       = "DC",
  "Monocyte" = "Monocyte",
  "NK_ILC"   = "NK/ILC",
  "NK/ILC"   = "NK/ILC",
  "gdT"      = "gdT",
  "MAIT"     = "MAIT",
  "Platelet" = "Platelet"
)

# Main populations shown in Figure 1.
# Treg is placed third, next to CD4 T and CD8 T.
cell_order_display <- c(
  "CD4 T",
  "CD8 T",
  "Treg",
  "B cell",
  "NK",
  "DC",
  "Monocyte"
)

ct_levels_switch <- c(
  "CD4_T",
  "CD8_T",
  "Treg",
  "B_cell",
  "NK",
  "DC",
  "Monocyte"
)

# Exclude minor/low-relevance populations consistently from B, C, D.
exclude_celltypes_raw <- c(
  "NK_ILC",
  "NK/ILC",
  "gdT",
  "Platelet",
  "MAIT",
  "pDC"
)

exclude_celltypes <- c(
  "NK/ILC",
  "gdT",
  "Platelet",
  "MAIT",
  "pDC"
)

# ============================================================
# TEXT SIZES
# ============================================================

base_size <- 14

# Main manuscript text size.
# This is the same size as the Panel B x-axis text.
axis_text_size <- 16

# Use the same size for all non-title text.
axis_title_size   <- axis_text_size
legend_text_size  <- axis_text_size
legend_title_size <- axis_text_size
subtitle_size     <- axis_text_size

# Titles and panel letters.
plot_title_size <- 16
panel_tag_size  <- 18

# ggplot2 geom_text() sizes are not in points.
# Convert 16 pt text to the corresponding geom_text size.
geom_text_size <- axis_text_size / ggplot2::.pt

panelA_bar_label_size <- geom_text_size
panelD_number_size    <- geom_text_size

# ============================================================
# PANEL B SUPPORT THRESHOLDS
# ============================================================

# A cell is counted as expressing the marker gene if the summed
# gene-level transcript count is at least 2.
# This is the version that preserves the 15-gene Panel B result.
cell_expr_threshold <- 2

# A gene-cell-type combination is supported if:
# 1) at least 10% of cells are expressing,
# 2) at least 10 cells are expressing,
# 3) the mean gene-level count among expressing cells is at least 2.
min_gene_expr_frac <- 0.10
min_expr_cells_abs <- 10L
min_mean_signal_per_expr_cell <- 2

# Keep genes for Panel B if at least one supported cell type has
# coding potential <= 0.5, meaning >=50% non-coding contribution.
gene_keep_threshold <- 0.50

# ============================================================
# PANEL C STRICT JS FILTERS
# ============================================================

min_pct_multitx_cells <- 10
min_expr_cells_strict <- 10L

# ============================================================
# COLORS
# ============================================================

col_text <- "grey15"

col_panelA <- "#4D4D4D"

col_grey_low  <- "#F2F2F2"
col_grey_high <- col_panelA

col_unsupported <- "#BDBDBD"

# Coding-potential gradient
col_noncoding <- "#E08A66"
col_mixed     <- "white"
col_coding    <- "#2D6A4F"

# ============================================================
# THEMES
# ============================================================

theme_fig <- ggplot2::theme_classic(base_size = base_size) +
  ggplot2::theme(
    text = ggplot2::element_text(color = col_text),
    axis.text = ggplot2::element_text(color = col_text, size = axis_text_size),
    axis.title = ggplot2::element_text(color = col_text, size = axis_title_size),
    plot.title = ggplot2::element_text(face = "bold", size = plot_title_size),
    plot.subtitle = ggplot2::element_text(size = subtitle_size, color = "grey35"),
    plot.tag = ggplot2::element_text(face = "bold", size = panel_tag_size),
    legend.text = ggplot2::element_text(size = legend_text_size),
    legend.title = ggplot2::element_text(size = legend_title_size),
    panel.grid = ggplot2::element_blank()
  )

theme_heat <- ggplot2::theme_classic(base_size = base_size) +
  ggplot2::theme(
    text = ggplot2::element_text(color = col_text),
    axis.text = ggplot2::element_text(color = col_text, size = axis_text_size),
    axis.title = ggplot2::element_text(color = col_text, size = axis_title_size),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = plot_title_size),
    plot.subtitle = ggplot2::element_text(size = subtitle_size, color = "grey35"),
    plot.caption = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(size = legend_title_size),
    legend.text = ggplot2::element_text(size = legend_text_size)
  )

# ============================================================
# HELPERS
# ============================================================

normalize_enst <- function(x) {
  stringr::str_extract(as.character(x), "ENST\\d+")
}

collapse_duplicate_rows <- function(mat) {
  ids <- rownames(mat)
  if (!anyDuplicated(ids)) return(mat)

  idx_list <- split(seq_len(nrow(mat)), ids)

  out_list <- lapply(names(idx_list), function(id) {
    ii <- idx_list[[id]]

    if (length(ii) == 1) {
      x <- mat[ii, , drop = FALSE]
      rownames(x) <- id
      return(x)
    } else {
      v <- Matrix::colSums(mat[ii, , drop = FALSE])
      x <- Matrix::Matrix(v, nrow = 1, sparse = TRUE)
      rownames(x) <- id
      colnames(x) <- colnames(mat)
      return(x)
    }
  })

  do.call(rbind, out_list)
}

get_assay_layer_safe <- function(obj, assay, layer) {
  out <- tryCatch(
    SeuratObject::LayerData(obj, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)

  out <- tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)

  out <- tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, slot = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)

  stop("Could not retrieve assay/layer: ", assay, " / ", layer)
}

make_corr_dist <- function(mat, margin = c("rows", "cols"), method = "pearson") {
  margin <- match.arg(margin)
  x <- if (margin == "rows") mat else t(mat)

  if (nrow(x) <= 1) {
    return(stats::dist(matrix(0, nrow = nrow(x), ncol = 1)))
  }

  cor_mat <- suppressWarnings(
    stats::cor(t(x), use = "pairwise.complete.obs", method = method)
  )

  cor_mat[is.na(cor_mat)] <- 0
  stats::as.dist(1 - cor_mat)
}

make_marker_labels <- function(label_df) {
  out <- label_df %>%
    dplyr::mutate(
      cd_name = stringr::str_trim(as.character(cd_name)),
      cd_name = dplyr::na_if(cd_name, ""),
      marker_label = dplyr::coalesce(cd_name, gene_symbol)
    )

  dup_idx <- duplicated(out$marker_label) | duplicated(out$marker_label, fromLast = TRUE)

  out$marker_label[dup_idx] <- ifelse(
    !is.na(out$cd_name[dup_idx]),
    paste0(out$cd_name[dup_idx], " (", out$gene_symbol[dup_idx], ")"),
    out$gene_symbol[dup_idx]
  )

  out
}

pretty_ct <- function(x) {
  out <- gsub("_", " ", x)
  out <- gsub("B cell", "B cell", out)
  out <- gsub("NK ILC", "NK/ILC", out)
  out
}

save_panel <- function(plot, filename_base, width, height, dpi = 600) {
  pdf_file <- file.path(fig_dir, paste0(filename_base, ".pdf"))
  png_file <- file.path(fig_dir, paste0(filename_base, ".png"))

  ggplot2::ggsave(
    pdf_file,
    plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )

  ggplot2::ggsave(
    png_file,
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )

  message("Saved: ", pdf_file)
  message("Saved: ", png_file)
}

# ============================================================
# LOAD INPUTS
# ============================================================

message("Loading transcript universe...")
tx_universe_pbmc <- readr::read_tsv(fig1_universe_file, show_col_types = FALSE)

message("Loading Seurat object...")
obj <- readRDS(seurat_file)

message("Loading CD marker metadata...")
cd_raw <- readr::read_csv(cd_marker_master_file, show_col_types = FALSE)

cd_gene_metadata <- cd_raw %>%
  dplyr::transmute(
    gene_symbol = as.character(.data$`Gene Symbol`),
    cd_name = as.character(.data$`CD Marker`)
  ) %>%
  dplyr::distinct()

# ============================================================
# PANEL A — TRANSCRIPT COMPLEXITY
# ============================================================

cap <- 50

gene_counts <- tx_universe_pbmc %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(
    transcripts_all = dplyr::n_distinct(ensembl_transcript_id),
    .groups = "drop"
  )

dist_tbl <- gene_counts %>%
  dplyr::mutate(
    tx_cap = dplyr::if_else(transcripts_all >= cap, cap, transcripts_all)
  ) %>%
  dplyr::count(tx_cap, name = "n_genes") %>%
  dplyr::arrange(tx_cap) %>%
  dplyr::mutate(
    xlab = dplyr::if_else(tx_cap == cap, paste0(cap, "+"), as.character(tx_cap))
  )

pA <- ggplot2::ggplot(
  dist_tbl,
  ggplot2::aes(x = factor(xlab, levels = xlab), y = n_genes)
) +
  ggplot2::geom_col(
    width = 0.82,
    fill = col_panelA,
    color = col_panelA,
    linewidth = 0.25
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_genes),
    vjust = -0.25,
    size = panelA_bar_label_size,
    color = col_text
  ) +
  ggplot2::scale_x_discrete(
    expand = ggplot2::expansion(add = 0.25)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.10))
  ) +
  ggplot2::labs(
    x = "PBMC-retained transcripts per gene",
    y = "CD marker genes",
    title = "Transcript complexity per CD marker gene"
  ) +
  theme_fig +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = axis_text_size),
    axis.text.y = ggplot2::element_text(size = axis_text_size),
    axis.title.x = ggplot2::element_text(size = axis_text_size),
    axis.title.y = ggplot2::element_text(size = axis_text_size),
    plot.margin = ggplot2::margin(5.5, 2, 5.5, 5.5)
  )

# ============================================================
# PANEL B — CELL-TYPE-SPECIFIC CODING POTENTIAL
# ============================================================

message("Preparing panel B...")

stopifnot(celltype_var %in% colnames(obj@meta.data))
stopifnot(isoform_assay %in% names(obj@assays))

mat <- get_assay_layer_safe(obj, isoform_assay, signal_layer)

rownames(mat) <- normalize_enst(rownames(mat))
mat <- mat[!is.na(rownames(mat)), , drop = FALSE]
mat <- collapse_duplicate_rows(mat)

celltype <- obj@meta.data[[celltype_var]]
names(celltype) <- rownames(obj@meta.data)

common_cells <- intersect(colnames(mat), names(celltype))
mat <- mat[, common_cells, drop = FALSE]
celltype <- celltype[common_cells]

keep_cells <- !is.na(celltype) &
  celltype != "" &
  !(celltype %in% exclude_celltypes_raw)

mat <- mat[, keep_cells, drop = FALSE]
celltype <- as.character(celltype[keep_cells])

celltype <- dplyr::recode(
  celltype,
  !!!celltype_map,
  .default = celltype
)

keep_cells2 <- !is.na(celltype) &
  celltype != "" &
  !(celltype %in% exclude_celltypes)

mat <- mat[, keep_cells2, drop = FALSE]
celltype <- celltype[keep_cells2]
names(celltype) <- colnames(mat)

tx_annot <- tx_universe_pbmc %>%
  dplyr::transmute(
    tx_id = normalize_enst(ensembl_transcript_id),
    gene_symbol = as.character(gene_symbol),
    transcript_biotype = as.character(transcript_biotype),
    coding_group = dplyr::case_when(
      transcript_biotype %in% c(
        "protein_coding",
        "protein_coding_CDS_not_defined",
        "protein_coding_LoF"
      ) ~ "Protein-coding-capable",
      TRUE ~ "Non-protein-coding"
    )
  ) %>%
  dplyr::filter(!is.na(tx_id), tx_id != "") %>%
  dplyr::distinct()

common_tx <- intersect(rownames(mat), tx_annot$tx_id)
mat <- mat[common_tx, , drop = FALSE]
tx_annot <- tx_annot %>%
  dplyr::filter(tx_id %in% common_tx)

agg_mat <- sapply(split(seq_along(celltype), celltype), function(ii) {
  Matrix::rowSums(mat[, ii, drop = FALSE])
})

if (is.null(dim(agg_mat))) {
  agg_mat <- matrix(
    agg_mat,
    ncol = 1,
    dimnames = list(names(agg_mat), unique(celltype))
  )
}

agg_long <- as.data.frame(agg_mat) %>%
  tibble::rownames_to_column("tx_id") %>%
  tidyr::pivot_longer(
    cols = -tx_id,
    names_to = "cell_type",
    values_to = "signal"
  )

tx_gene_map <- tx_annot %>%
  dplyr::distinct(tx_id, gene_symbol)

common_tx_support <- intersect(rownames(mat), tx_gene_map$tx_id)

tx_gene_map <- tx_gene_map %>%
  dplyr::filter(tx_id %in% common_tx_support)

tx_gene_map <- tx_gene_map[match(common_tx_support, tx_gene_map$tx_id), , drop = FALSE]
mat_support <- mat[common_tx_support, , drop = FALSE]

gene_levels <- unique(tx_gene_map$gene_symbol)

G <- Matrix::sparseMatrix(
  i = match(tx_gene_map$gene_symbol, gene_levels),
  j = seq_along(common_tx_support),
  x = 1,
  dims = c(length(gene_levels), length(common_tx_support)),
  dimnames = list(gene_levels, common_tx_support)
)

gene_cell_mat <- G %*% mat_support

support_tbl <- dplyr::bind_rows(lapply(split(seq_along(celltype), celltype), function(ii) {
  ct <- unique(celltype[ii])
  sub_mat <- gene_cell_mat[, ii, drop = FALSE]

  n_cells_expr_vec <- as.numeric(Matrix::rowSums(sub_mat >= cell_expr_threshold))
  gene_total_signal_vec <- as.numeric(Matrix::rowSums(sub_mat))

  tibble::tibble(
    gene_symbol = rownames(sub_mat),
    cell_type = ct,
    n_cells_total = length(ii),
    n_cells_expr = n_cells_expr_vec,
    frac_cells_expr = n_cells_expr_vec / length(ii),
    gene_total_signal = gene_total_signal_vec,
    mean_signal_per_cell = gene_total_signal_vec / length(ii),
    mean_signal_per_expr_cell = dplyr::if_else(
      n_cells_expr_vec > 0,
      gene_total_signal_vec / n_cells_expr_vec,
      NA_real_
    )
  )
}))

plot_tbl_B_all <- agg_long %>%
  dplyr::inner_join(
    tx_annot %>% dplyr::select(tx_id, gene_symbol, coding_group),
    by = "tx_id"
  ) %>%
  dplyr::group_by(gene_symbol, cell_type, coding_group) %>%
  dplyr::summarise(signal = sum(signal, na.rm = TRUE), .groups = "drop") %>%
  tidyr::complete(
    gene_symbol,
    cell_type,
    coding_group = c("Protein-coding-capable", "Non-protein-coding"),
    fill = list(signal = 0)
  ) %>%
  dplyr::group_by(gene_symbol, cell_type) %>%
  dplyr::summarise(
    total_signal = sum(signal),
    coding_signal = sum(signal[coding_group == "Protein-coding-capable"]),
    coding_potential_raw = dplyr::if_else(
      total_signal > 0,
      coding_signal / total_signal,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(support_tbl, by = c("gene_symbol", "cell_type")) %>%
  dplyr::mutate(
    pass_expression_frac = frac_cells_expr >= min_gene_expr_frac,
    pass_expression_abs = n_cells_expr >= min_expr_cells_abs,
    pass_signal_strength = mean_signal_per_expr_cell >= min_mean_signal_per_expr_cell,
    pass_support = pass_expression_frac &
      pass_expression_abs &
      pass_signal_strength,
    coding_potential = dplyr::if_else(
      pass_support,
      coding_potential_raw,
      NA_real_
    )
  )

readr::write_csv(
  plot_tbl_B_all,
  file.path(fig_dir, "Figure1B_coding_potential_support_full_unfiltered.csv")
)

readr::write_csv(
  plot_tbl_B_all,
  file.path(fig_dir, "Figure1B_coding_potential_support_10pct_cells_count2.csv")
)

# Gene inclusion for Panel B.
# This is computed from the full unfiltered table to avoid stale/filtered-state errors.
gene_keep_tbl_B <- plot_tbl_B_all %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(
    keep_gene = any(
      pass_support & coding_potential_raw <= gene_keep_threshold,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

message("Number of Panel B genes shown: ", sum(gene_keep_tbl_B$keep_gene))

panelB_kept_genes <- gene_keep_tbl_B %>%
  dplyr::filter(keep_gene) %>%
  dplyr::left_join(cd_gene_metadata, by = "gene_symbol") %>%
  dplyr::mutate(
    marker_label = dplyr::coalesce(cd_name, gene_symbol)
  ) %>%
  dplyr::arrange(marker_label) %>%
  dplyr::select(marker_label, gene_symbol)

print(panelB_kept_genes, n = Inf)

readr::write_csv(
  panelB_kept_genes,
  file.path(fig_dir, "Figure1B_kept_genes_METHODS_filter.csv")
)

# Diagnostic comparison to the older total-signal support rule.
old_vs_methods_gene_keep <- plot_tbl_B_all %>%
  dplyr::mutate(
    old_min_cells_required = pmax(10L, ceiling(0.02 * n_cells_total)),
    old_pass_support = total_signal >= 100 &
      n_cells_expr >= old_min_cells_required,
    methods_pass_support = pass_support,
    old_keep_tile = old_pass_support &
      coding_potential_raw <= gene_keep_threshold,
    methods_keep_tile = methods_pass_support &
      coding_potential_raw <= gene_keep_threshold
  ) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(
    old_keep = any(old_keep_tile, na.rm = TRUE),
    methods_keep = any(methods_keep_tile, na.rm = TRUE),
    old_supporting_celltypes = paste(
      sort(unique(cell_type[old_keep_tile])),
      collapse = ", "
    ),
    methods_supporting_celltypes = paste(
      sort(unique(cell_type[methods_keep_tile])),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(cd_gene_metadata, by = "gene_symbol") %>%
  dplyr::mutate(
    marker_label = dplyr::coalesce(cd_name, gene_symbol),
    status = dplyr::case_when(
      old_keep & methods_keep ~ "kept_both",
      old_keep & !methods_keep ~ "old_only",
      !old_keep & methods_keep ~ "methods_only",
      TRUE ~ "kept_neither"
    )
  )

message("Old total-signal rule versus Methods-style rule:")
print(
  old_vs_methods_gene_keep %>%
    dplyr::count(status),
  n = Inf
)

readr::write_csv(
  old_vs_methods_gene_keep,
  file.path(fig_dir, "Figure1B_old_total100_vs_methods_support_gene_keep_FULL.csv")
)

# Filter for plotting only after the full support and diagnostic tables are saved.
plot_tbl_B <- plot_tbl_B_all %>%
  dplyr::inner_join(
    gene_keep_tbl_B %>%
      dplyr::filter(keep_gene) %>%
      dplyr::select(gene_symbol),
    by = "gene_symbol"
  )

label_tbl_B <- plot_tbl_B %>%
  dplyr::distinct(gene_symbol) %>%
  dplyr::left_join(cd_gene_metadata, by = "gene_symbol") %>%
  make_marker_labels()

plot_tbl_B <- plot_tbl_B %>%
  dplyr::left_join(
    label_tbl_B %>% dplyr::select(gene_symbol, cd_name, marker_label),
    by = "gene_symbol"
  )

ord_B <- plot_tbl_B %>%
  dplyr::select(gene_symbol, cell_type, coding_potential) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    names_from = cell_type,
    values_from = coding_potential
  ) %>%
  tibble::column_to_rownames("gene_symbol")

ord_mat_B <- as.matrix(ord_B)

gene_dist_B <- make_corr_dist(ord_mat_B, margin = "rows", method = "pearson")
gene_order_B <- rownames(ord_mat_B)[stats::hclust(gene_dist_B, method = "average")$order]

marker_order_B <- label_tbl_B %>%
  dplyr::filter(gene_symbol %in% gene_order_B) %>%
  dplyr::mutate(gene_rank = match(gene_symbol, gene_order_B)) %>%
  dplyr::arrange(gene_rank) %>%
  dplyr::pull(marker_label)

plot_tbl_B <- plot_tbl_B %>%
  dplyr::mutate(
    marker_label = factor(marker_label, levels = marker_order_B),
    cell_type = factor(cell_type, levels = rev(cell_order_display))
  )

unsupported_B <- plot_tbl_B %>% dplyr::filter(!pass_support)
supported_B   <- plot_tbl_B %>% dplyr::filter(pass_support)

pB <- ggplot2::ggplot() +
  ggplot2::geom_tile(
    data = unsupported_B,
    ggplot2::aes(x = marker_label, y = cell_type),
    fill = col_unsupported,
    colour = NA
  ) +
  ggplot2::geom_tile(
    data = supported_B,
    ggplot2::aes(x = marker_label, y = cell_type, fill = coding_potential),
    colour = NA
  ) +
  ggplot2::scale_fill_gradient2(
    low = col_noncoding,
    mid = col_mixed,
    high = col_coding,
    midpoint = 0.5,
    limits = c(0, 1),
    oob = scales::squish,
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0\nnon-coding", "0.25", "0.5\nmixed", "0.75", "1\ncoding"),
    name = "Coding\npotential\n"
  ) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0)) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0)) +
  ggplot2::labs(
    x = "",
    y = "Cell type",
    title = "Cell-type-specific coding potential"
  ) +
  theme_heat +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 0.5,
      size = axis_text_size
    ),
    axis.text.y = ggplot2::element_text(size = axis_text_size),
    axis.title.y = ggplot2::element_text(size = axis_text_size),
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank(),
    legend.position = "right"
  )

# ============================================================
# PANEL C — STRICT WITHIN-POPULATION JS BOXPLOT
# ============================================================

message("Preparing panel C...")

if (!file.exists(js_file)) {
  stop("Missing JS file: ", js_file, "\nRun the within-population JS calculation first.")
}

if (!file.exists(diag_file)) {
  stop("Missing diagnostic file: ", diag_file, "\nRun the JS sparsity diagnostic first.")
}

js_tbl <- readr::read_csv(js_file, show_col_types = FALSE)
diag_tbl <- readr::read_csv(diag_file, show_col_types = FALSE)

js_tbl <- js_tbl %>%
  dplyr::mutate(
    cell_type = dplyr::recode(
      as.character(cell_type),
      !!!celltype_map,
      .default = as.character(cell_type)
    )
  )

diag_tbl <- diag_tbl %>%
  dplyr::mutate(
    cell_type = dplyr::recode(
      as.character(cell_type),
      !!!celltype_map,
      .default = as.character(cell_type)
    )
  )

diag_tbl <- diag_tbl %>%
  dplyr::select(
    gene_symbol,
    cell_type,
    dplyr::any_of(c(
      "n_expr_cells",
      "median_detected_pc_tx_per_cell",
      "mean_detected_pc_tx_per_cell",
      "pct_expr_cells_with_one_pc_tx",
      "pct_expr_cells_with_two_or_more_pc_tx"
    ))
  ) %>%
  dplyr::distinct(gene_symbol, cell_type, .keep_all = TRUE)

label_tbl_C <- js_tbl %>%
  dplyr::distinct(gene_symbol) %>%
  dplyr::left_join(cd_gene_metadata, by = "gene_symbol") %>%
  make_marker_labels()

heat_base_C <- js_tbl %>%
  dplyr::left_join(diag_tbl, by = c("gene_symbol", "cell_type")) %>%
  dplyr::left_join(
    label_tbl_C %>% dplyr::select(gene_symbol, cd_name, marker_label),
    by = "gene_symbol"
  )

if (!"n_expr_cells" %in% colnames(heat_base_C)) {
  heat_base_C$n_expr_cells <- heat_base_C$n_cells_pc_expr
} else if ("n_cells_pc_expr" %in% colnames(heat_base_C)) {
  heat_base_C$n_expr_cells <- dplyr::coalesce(
    heat_base_C$n_expr_cells,
    heat_base_C$n_cells_pc_expr
  )
}

heat_base_C <- heat_base_C %>%
  dplyr::mutate(
    strict_pass = pass_support &
      !is.na(median_js_distance) &
      !is.na(pct_expr_cells_with_two_or_more_pc_tx) &
      pct_expr_cells_with_two_or_more_pc_tx >= min_pct_multitx_cells &
      n_expr_cells >= min_expr_cells_strict,

    support_class = dplyr::case_when(
      strict_pass ~ "Strict-supported",
      !pass_support ~ "Insufficient expression support",
      is.na(pct_expr_cells_with_two_or_more_pc_tx) ~ "No multi-transcript diagnostic",
      pct_expr_cells_with_two_or_more_pc_tx < min_pct_multitx_cells ~ "Fails multi-transcript support",
      n_expr_cells < min_expr_cells_strict ~ "Too few expressing cells",
      TRUE ~ "Other filtered"
    )
  )

strict_panel_summary <- heat_base_C %>%
  dplyr::filter(
    strict_pass,
    cell_type %in% cell_order_display
  ) %>%
  dplyr::summarise(
    n_strict_gene_celltypes = dplyr::n(),
    n_strict_genes = dplyr::n_distinct(gene_symbol),
    n_celltypes = dplyr::n_distinct(cell_type),
    median_js = median(median_js_distance, na.rm = TRUE),
    mean_js = mean(median_js_distance, na.rm = TRUE),
    pct_gene_celltypes_js_ge_0_50 =
      100 * mean(median_js_distance >= 0.50, na.rm = TRUE),
    pct_gene_celltypes_js_ge_0_75 =
      100 * mean(median_js_distance >= 0.75, na.rm = TRUE),
    median_pct_multitx_cells =
      median(pct_expr_cells_with_two_or_more_pc_tx, na.rm = TRUE)
  )

print(strict_panel_summary)

readr::write_csv(
  strict_panel_summary,
  file.path(fig_dir, "Figure1C_strict_JS_boxplot_summary.csv")
)

readr::write_csv(
  heat_base_C,
  file.path(fig_dir, "Figure1C_strict_JS_boxplot_full_table.csv")
)

strict_js_plot_tbl <- heat_base_C %>%
  dplyr::filter(
    strict_pass,
    cell_type %in% cell_order_display
  ) %>%
  dplyr::mutate(
    cell_type = factor(cell_type, levels = cell_order_display)
  )

pC <- ggplot2::ggplot(
  strict_js_plot_tbl,
  ggplot2::aes(x = cell_type, y = median_js_distance)
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "grey92",
    colour = "grey25",
    linewidth = 0.35
  ) +
  ggplot2::geom_jitter(
    width = 0.16,
    height = 0,
    size = 0.75,
    alpha = 0.45,
    colour = "grey35"
  ) +
  ggplot2::geom_hline(
    yintercept = c(0.5, 0.75),
    linetype = "dashed",
    linewidth = 0.3,
    colour = "grey55"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = ggplot2::expansion(mult = c(0.02, 0.04))
  ) +
  ggplot2::labs(
    x = "Cell type",
    y = "Median pairwise JS distance",
    title = "Within-population coding transcript heterogeneity"
  ) +
  theme_fig +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = axis_text_size
    ),
    axis.text.y = ggplot2::element_text(size = axis_text_size),
    axis.title.x = ggplot2::element_text(
      size = axis_text_size,
      margin = ggplot2::margin(t = 10)
    ),
    axis.title.y = ggplot2::element_text(size = axis_text_size),
    plot.caption = ggplot2::element_text(size = 8, color = "grey35", hjust = 0),
    legend.position = "none"
  )

# ============================================================
# PANEL D — PAIRWISE ISOFORM SWITCHING TRIANGLE
# ============================================================

message("Preparing panel D...")

gene_calls_pc <- readRDS(gene_calls_pc_file)
usable_ct_pc  <- readRDS(usable_ct_pc_file)

build_switch_triangle <- function(gene_calls, usable_ct, ct_levels) {
  ct_order <- intersect(ct_levels, usable_ct)

  pair_counts <- gene_calls %>%
    dplyr::filter(
      condition_1 %in% ct_order,
      condition_2 %in% ct_order
    ) %>%
    dplyr::mutate(
      i1 = match(condition_1, ct_order),
      i2 = match(condition_2, ct_order)
    ) %>%
    dplyr::filter(!is.na(i1), !is.na(i2), i1 != i2) %>%
    dplyr::mutate(
      ct_left  = ifelse(i1 < i2, condition_1, condition_2),
      ct_right = ifelse(i1 < i2, condition_2, condition_1)
    ) %>%
    dplyr::count(ct_left, ct_right, name = "n_switching_genes") %>%
    dplyr::mutate(
      x = factor(ct_left, levels = ct_order, labels = pretty_ct(ct_order)),
      y = factor(ct_right, levels = rev(ct_order), labels = pretty_ct(rev(ct_order))),
      label_col = ifelse(n_switching_genes >= 40, "white", "grey15")
    )

  ggplot2::ggplot(pair_counts, ggplot2::aes(x = x, y = y, fill = n_switching_genes)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_switching_genes, color = label_col),
      size = panelD_number_size
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient(
      low = col_grey_low,
      high = col_grey_high,
      name = "# genes with\nisoform\nswitch"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Isoform switching across populations"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        color = col_text,
        size = axis_text_size
      ),
      axis.text.y = ggplot2::element_text(
        color = col_text,
        size = axis_text_size
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = plot_title_size,
        color = col_text
      ),
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size)
    ) +
    ggplot2::coord_equal()
}

pD <- build_switch_triangle(
  gene_calls = gene_calls_pc,
  usable_ct = usable_ct_pc,
  ct_levels = ct_levels_switch
)

# ============================================================
# SAVE SEPARATE PANELS
# ============================================================

save_panel(
  pA,
  "Figure1A_transcript_complexity",
  width = 13.5,
  height = 3.3,
  dpi = 600
)

save_panel(
  pB,
  "Figure1B_coding_potential_heatmap",
  width = 11.5,
  height = 4.6,
  dpi = 600
)

save_panel(
  pC,
  "Figure1C_strict_JS_boxplot",
  width = 6.4,
  height = 4.6,
  dpi = 600
)

save_panel(
  pD,
  "Figure1D_pairwise_isoform_switching",
  width = 5.4,
  height = 5.0,
  dpi = 600
)

# ============================================================
# COMBINED FIGURE
# ============================================================

bottom_row <- pC | pD

pA_for_combined <- patchwork::free(pA, side = "r")

fig_all <- pA_for_combined / pB / bottom_row +
  patchwork::plot_layout(
    heights = c(1.00, 1.35, 1.48),
    widths = c(1.65, 0.85)
  ) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = panel_tag_size,
      color = col_text
    ),
    plot.margin = ggplot2::margin(5.5, 7, 5.5, 5.5)
  )

ggplot2::ggsave(
  outfile_combined_pdf,
  fig_all,
  width = 18.0,
  height = 13.4,
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  outfile_combined_png,
  fig_all,
  width = 18.0,
  height = 13.4,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

message("Saved combined PDF: ", outfile_combined_pdf)
message("Saved combined PNG: ", outfile_combined_png)
