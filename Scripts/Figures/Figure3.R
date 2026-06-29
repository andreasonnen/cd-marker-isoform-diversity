#!/usr/bin/env Rscript

# ============================================================
# Generate Figure 3: CD45 isoform usage, estimated molecular
# weight, marker expression, and functional associations.

# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(msigdbr)
  library(AUCell)
  library(scales)
  library(RColorBrewer)
})

set.seed(1)

# ============================================================
# 0) CONFIG
# ============================================================

input_rds <- Sys.getenv(
  "CD_SEURAT_RDS",
  unset = "../Data/combined_runs.rds"
)

output_dir <- Sys.getenv(
  "CD_FIGURE_DIR",
  unset = "../Figures"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

celltype_col <- "broad_celltype_final"
celltype_val <- "CD8_T"

isoform_assay <- "Isoform"
gene_assay <- "Gene"

isoform_tpm_layer <- "tpm"
gene_expr_layer <- "data"
gene_tpm_layer <- "tpm"

mw_col <- "CD45_MW_kDa_weighted"

target_cells_per_bin <- 30
min_total_cd45_tpm <- 0.5
seed_bins <- 1

clip_marker_z <- 3
clip_path_z <- 3

fig_width <- 8
fig_height <- 10
fig_dpi <- 600

outfile_pdf <- file.path(
  output_dir,
  "Figure3_CD45_isoform_story_domIsoBins_domThenMW.pdf"
)

outfile_png <- file.path(
  output_dir,
  "Figure3_CD45_isoform_story_domIsoBins_domThenMW.png"
)

iso_order <- c("RABC", "RAB", "RBC", "RA", "RB", "RO")

iso_labels <- c(
  RABC = "CD45RABC",
  RAB  = "CD45RAB",
  RBC  = "CD45RBC",
  RA   = "CD45RA",
  RB   = "CD45RB",
  RO   = "CD45RO"
)

marker_genes <- c(
  "CCR7", "SELL", "IL7R", "LTB", "LEF1", "TCF7",
  "NKG7", "PRF1", "GZMB", "CCL5", "GNLY"
)

pw_keep <- c(
  "GOBP_CYTOPLASMIC_TRANSLATION",
  "GOBP_RIBOSOME_BIOGENESIS",
  "GOBP_RNA_SPLICING",
  "GOBP_POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_PROLIFERATION",
  "GOBP_REGULATION_OF_T_CELL_ACTIVATION",
  "GOBP_LYMPHOCYTE_DIFFERENTIATION",
  "GOBP_REGULATION_OF_LAMELLIPODIUM_ASSEMBLY",
  "GOBP_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS",
  "GOBP_FC_RECEPTOR_MEDIATED_STIMULATORY_SIGNALING_PATHWAY",
  "GOBP_T_CELL_MEDIATED_CYTOTOXICITY"
)

# ============================================================
# 1) HELPERS
# ============================================================

strip_version <- function(x) {
  sub("\\.\\d+$", "", x)
}

get_layer_safe <- function(seu, assay, layer, slot_fallback = NULL) {
  if (!assay %in% names(seu@assays)) {
    return(NULL)
  }

  out <- tryCatch(
    SeuratObject::LayerData(seu, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)

  out <- tryCatch(
    SeuratObject::GetAssayData(seu, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)

  if (!is.null(slot_fallback)) {
    out <- tryCatch(
      SeuratObject::GetAssayData(seu, assay = assay, slot = slot_fallback),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }

  NULL
}

join_layer_if_possible <- function(obj, assay, layer) {
  if (!assay %in% names(obj@assays)) {
    return(obj)
  }

  obj[[assay]] <- tryCatch(
    {
      SeuratObject::JoinLayers(obj[[assay]], layers = layer)
    },
    error = function(e) {
      message("Skipping JoinLayers for ", assay, " / ", layer, ": ", e$message)
      obj[[assay]]
    }
  )

  obj
}

summarise_auc_by_bin <- function(auc_mat, bin_df) {
  stopifnot(all(c("cell", "bin") %in% colnames(bin_df)))

  common <- intersect(bin_df$cell, colnames(auc_mat))
  if (length(common) == 0) {
    stop("No overlap between bin_df cells and auc_mat columns.")
  }

  dfb <- bin_df %>%
    dplyr::filter(cell %in% common)

  bins <- sort(unique(dfb$bin))

  pw_mean <- sapply(bins, function(b) {
    cells_b <- dfb$cell[dfb$bin == b]
    rowMeans(auc_mat[, cells_b, drop = FALSE], na.rm = TRUE)
  })

  pw_mean <- as.matrix(pw_mean)
  colnames(pw_mean) <- as.character(bins)

  pw_mean
}

# ============================================================
# 2) CD45 MW + CANONICAL ISOFORM EXPRESSION
# ============================================================

compute_cd45_mw <- function(
    obj,
    assay = "Isoform",
    tpm_layer = "tpm",
    mw_col = "CD45_MW_kDa_weighted",
    iso_levels = c("RO", "RB", "RA", "RAB", "RBC", "RABC"),
    N_OCC = 0.85,
    O_OCC = 0.60,
    N_PER = 2.40,
    O_PER = 0.95
) {
  iso_mass_tbl <- tibble::tribble(
    ~isoform, ~poly_kDa, ~N_sites, ~O_sites,
    "RO",     131.130,   13,       32,
    "RB",     135.719,   12,       39,
    "RA",     138.034,   15,       41,
    "RAB",    142.623,   15,       60,
    "RBC",    140.582,   14,       50,
    "RABC",   147.486,   17,       75
  ) %>%
    dplyr::mutate(
      N_added_kDa = N_sites * N_OCC * N_PER,
      O_added_kDa = O_sites * O_OCC * O_PER,
      glycan_mass_kDa = N_added_kDa + O_added_kDa,
      apparent_kDa = poly_kDa + glycan_mass_kDa
    ) %>%
    dplyr::filter(isoform %in% iso_levels)

  tpm_mat <- get_layer_safe(
    seu = obj,
    assay = assay,
    layer = tpm_layer,
    slot_fallback = tpm_layer
  )

  if (is.null(tpm_mat)) {
    stop("Could not read TPM from assay = ", assay, ", layer = ", tpm_layer)
  }

  stopifnot(ncol(tpm_mat) == ncol(obj))

  # Final transcript-to-CD45-isoform mapping used for Figure 3.
  # Transcript versions are stripped before matching, so this is robust
  # to version suffix differences in the assay rownames.
  cd45_tx_to_isoform <- tibble::tribble(
    ~feature,             ~isoform,
    "ENST00000442510.8",  "RABC",
    "ENST00000697631.1",  "RA",
    "ENST00000529828.5",  "RAB",
    "ENST00000970625.1",  "RB",
    "ENST00000367367.8",  "RBC",
    "ENST00000348564.12", "RO",
    "ENST00000367379.6",  "RO",
    "ENST00000697632.1",  "RO",
    "ENST00000908298.1",  "RO",
    "ENST00000908299.1",  "RO",
    "ENST00000970623.1",  "RO",
    "ENST00000970624.1",  "RO",
    "ENST00000970626.1",  "RO"
  )

  feature_df <- tibble::tibble(feature = rownames(tpm_mat)) %>%
    dplyr::mutate(feature_novers = strip_version(feature))

  cd45_tx_to_isoform <- cd45_tx_to_isoform %>%
    dplyr::mutate(feature_novers = strip_version(feature))

  map_in_assay <- cd45_tx_to_isoform %>%
    dplyr::inner_join(feature_df, by = "feature_novers", suffix = c(".map", ".assay")) %>%
    dplyr::transmute(
      feature = feature.assay,
      isoform = isoform
    ) %>%
    dplyr::filter(isoform %in% iso_levels) %>%
    dplyr::distinct()

  if (nrow(map_in_assay) == 0) {
    stop("No canonical CD45 transcripts found in the Isoform TPM assay.")
  }

  message("CD45 transcript-to-isoform mapping found in assay:")
  print(map_in_assay, n = Inf)

  feat2iso <- stats::setNames(map_in_assay$isoform, map_in_assay$feature)

  iso_expr <- sapply(iso_levels, function(iso) {
    feats <- names(feat2iso)[feat2iso == iso]
    feats <- intersect(feats, rownames(tpm_mat))

    if (length(feats) == 0) {
      return(rep(0, ncol(tpm_mat)))
    }

    Matrix::colSums(tpm_mat[feats, , drop = FALSE])
  })

  iso_expr <- as.matrix(iso_expr)
  rownames(iso_expr) <- colnames(tpm_mat)
  stopifnot(all(colnames(iso_expr) == iso_levels))

  mass_vec <- iso_mass_tbl$apparent_kDa[
    match(colnames(iso_expr), iso_mass_tbl$isoform)
  ]
  names(mass_vec) <- colnames(iso_expr)

  if (any(!is.finite(mass_vec))) {
    stop("mass_vec contains NA/Inf.")
  }

  w_sum <- rowSums(iso_expr)

  mw_vec <- as.numeric((iso_expr %*% mass_vec) / pmax(w_sum, 1e-12))
  mw_vec[w_sum == 0] <- NA_real_
  names(mw_vec) <- rownames(iso_expr)

  obj <- SeuratObject::AddMetaData(
    object = obj,
    metadata = mw_vec,
    col.name = mw_col
  )

  list(
    obj = obj,
    iso_expr = iso_expr,
    iso_mass_tbl = iso_mass_tbl,
    map_in_assay = map_in_assay
  )
}

# ============================================================
# 3) DOMINANT-ISOFORM BINNING
# ============================================================

make_bins_by_dominant_isoform <- function(
    iso_expr_cells_by_iso,
    mw_vec = NULL,
    target_cells_per_bin = 30,
    dom_order = c("RABC", "RAB", "RBC", "RA", "RB", "RO"),
    min_total_cd45_tpm = 0.5,
    seed = 1
) {
  stopifnot(is.matrix(iso_expr_cells_by_iso) || inherits(iso_expr_cells_by_iso, "Matrix"))

  dom_order <- intersect(dom_order, colnames(iso_expr_cells_by_iso))
  if (length(dom_order) == 0) {
    stop("No requested isoforms found in isoform expression matrix.")
  }

  X <- iso_expr_cells_by_iso[, dom_order, drop = FALSE]

  keep <- rowSums(X) >= min_total_cd45_tpm
  X <- X[keep, , drop = FALSE]

  if (nrow(X) == 0) {
    stop("No cells pass min_total_cd45_tpm.")
  }

  frac <- X / pmax(rowSums(X), 1e-12)
  frac <- as.matrix(frac)

  dom_idx <- max.col(frac, ties.method = "first")
  dom_iso <- dom_order[dom_idx]
  dom_frac <- apply(frac, 1, max)

  df <- tibble::tibble(
    cell = rownames(frac),
    dominant_isoform = factor(dom_iso, levels = dom_order),
    dom_frac = dom_frac
  )

  if (!is.null(mw_vec)) {
    df$MW <- as.numeric(mw_vec[df$cell])
  } else {
    df$MW <- NA_real_
  }

  set.seed(seed)

  df <- df %>%
    dplyr::mutate(tie_rand = runif(dplyr::n()))

  out_list <- list()
  bin_counter <- 1L

  for (iso in dom_order) {
    sub <- df %>%
      dplyr::filter(dominant_isoform == iso) %>%
      dplyr::arrange(dplyr::desc(dom_frac), dplyr::desc(MW), tie_rand)

    if (nrow(sub) == 0) next

    n_bins_iso <- ceiling(nrow(sub) / target_cells_per_bin)

    base <- floor(nrow(sub) / n_bins_iso)
    rem <- nrow(sub) %% n_bins_iso
    sizes <- c(rep(base + 1, rem), rep(base, n_bins_iso - rem))

    sub$bin <- rep(seq(bin_counter, length.out = n_bins_iso), times = sizes)
    bin_counter <- max(sub$bin) + 1L

    out_list[[iso]] <- sub
  }

  bin_df <- dplyr::bind_rows(out_list) %>%
    dplyr::select(cell, dominant_isoform, dom_frac, MW, bin)

  list(
    bin_df = bin_df,
    iso_frac_per_cell = frac
  )
}

# ============================================================
# 4) MARKER HEATMAP HELPER
# ============================================================

make_marker_bin_heatmap <- function(
    obj_sub,
    bin_df,
    markers,
    assay_priority = c("Gene", "RNA"),
    layer = "tpm",
    which = c("mean", "pct"),
    pct_threshold = 0,
    bin_levels = NULL,
    clip_z = 3
) {
  which <- match.arg(which)

  gene_mat <- NULL

  for (a in assay_priority) {
    gene_mat <- get_layer_safe(obj_sub, assay = a, layer = layer, slot_fallback = layer)
    if (!is.null(gene_mat)) break
  }

  if (is.null(gene_mat)) {
    stop("Could not read expression matrix from requested Gene/RNA assay layers.")
  }

  markers_present <- intersect(markers, rownames(gene_mat))

  if (length(markers_present) < 3) {
    stop("Too few marker genes found. Found: ", paste(markers_present, collapse = ", "))
  }

  common_cells <- intersect(bin_df$cell, colnames(gene_mat))

  bin_df_use <- bin_df %>%
    dplyr::filter(cell %in% common_cells)

  bin_vec <- bin_df_use$bin
  names(bin_vec) <- bin_df_use$cell

  bins <- sort(unique(bin_vec))

  marker_bin <- sapply(bins, function(b) {
    cells_b <- names(bin_vec)[bin_vec == b]
    x <- gene_mat[markers_present, cells_b, drop = FALSE]

    if (which == "mean") {
      Matrix::rowMeans(x)
    } else {
      Matrix::rowMeans(x > pct_threshold)
    }
  })

  marker_bin <- as.matrix(marker_bin)
  rownames(marker_bin) <- markers_present
  colnames(marker_bin) <- as.character(bins)

  z <- t(scale(t(marker_bin)))
  z[!is.finite(z)] <- 0

  gene_clust <- hclust(dist(z))
  gene_order <- rownames(z)[gene_clust$order]

  df <- as.data.frame(z) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(
      cols = -gene,
      names_to = "bin",
      values_to = "z"
    ) %>%
    dplyr::mutate(
      gene = factor(gene, levels = gene_order),
      bin = as.integer(bin)
    )

  if (!is.null(bin_levels)) {
    df$bin <- factor(df$bin, levels = bin_levels)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = bin, y = gene, fill = z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
      limits = c(-clip_z, clip_z),
      values = scales::rescale(seq(-clip_z, clip_z, length.out = 100)),
      breaks = c(-clip_z, clip_z),
      oob = scales::squish,
      name = "Mean (z)"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Gene expression"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(hjust = 0, face = "bold", size = 10)
    )

  list(
    plot = p,
    marker_bin = marker_bin,
    z = z,
    markers_present = markers_present,
    bin_df_use = bin_df_use
  )
}

# ============================================================
# 5) PLOT HELPERS
# ============================================================

plot_mw <- function(
    mw_sum,
    bin_levels,
    title = "Median estimated MW (kDa)",
    hide_x = TRUE
) {
  df <- mw_sum %>%
    dplyr::mutate(bin_f = factor(bin, levels = bin_levels)) %>%
    dplyr::arrange(bin_f)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = bin_f, y = medMW, group = 1)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title.position = "panel",
      plot.title = ggplot2::element_text(hjust = 0, face = "bold", size = 10),
      axis.text = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(t = 2, r = 5.5, b = 5.5, l = 5.5, unit = "pt")
    )

  if (hide_x) {
    p <- p +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  p
}

make_story_plot <- function(
    bin_levels,
    iso_frac,
    mw_sum,
    pathway_z_mat,
    obj,
    bin_df,
    marker_genes,
    clip_marker_z = 3,
    clip_path_z = 3,
    iso_order = c("RABC", "RAB", "RBC", "RA", "RB", "RO"),
    iso_labels,
    pathway_order
) {
  yl_or_rd_w <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100)
  rd_bu_rev_w <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)

  legend_barwidth <- grid::unit(0.4, "cm")
  legend_barheight <- grid::unit(2.5, "cm")

  legend_theme <- ggplot2::theme(
    legend.position = "right",
    legend.justification = c(0, 0.5),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 0, 2, unit = "pt")
  )

  tight_margin_top <- ggplot2::theme(
    plot.margin = ggplot2::margin(t = 2, r = 5.5, b = 1, l = 5.5, unit = "pt")
  )

  tight_margin_mid <- ggplot2::theme(
    plot.margin = ggplot2::margin(t = 1, r = 5.5, b = 1, l = 5.5, unit = "pt")
  )

  tight_margin_bot <- ggplot2::theme(
    plot.margin = ggplot2::margin(t = 1, r = 5.5, b = 5.5, l = 5.5, unit = "pt")
  )

  # Panel A: isoform usage
  df_iso <- as.data.frame(iso_frac) %>%
    tibble::rownames_to_column("isoform") %>%
    tidyr::pivot_longer(
      cols = -isoform,
      names_to = "bin",
      values_to = "frac"
    ) %>%
    dplyr::mutate(
      bin = factor(as.integer(bin), levels = bin_levels),
      isoform = factor(isoform, levels = rev(iso_order))
    )

  p_iso <- ggplot2::ggplot(df_iso, ggplot2::aes(x = bin, y = isoform, fill = frac)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = yl_or_rd_w,
      limits = c(0, 1),
      oob = scales::squish,
      name = "Fraction"
    ) +
    ggplot2::scale_y_discrete(labels = iso_labels) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        barwidth = legend_barwidth,
        barheight = legend_barheight
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "CD45 isoform usage"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(hjust = 0, face = "bold", size = 10)
    ) +
    legend_theme +
    tight_margin_top

  # Panel B: marker heatmap
  obj_cells <- intersect(bin_df$cell, colnames(obj))
  obj_mark <- subset(obj, cells = obj_cells)

  mk <- make_marker_bin_heatmap(
    obj_sub = obj_mark,
    bin_df = bin_df,
    markers = marker_genes,
    assay_priority = c("Gene", "RNA"),
    layer = gene_tpm_layer,
    which = "mean",
    bin_levels = bin_levels,
    clip_z = clip_marker_z
  )

  p_markers <- mk$plot +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        barwidth = legend_barwidth,
        barheight = legend_barheight
      )
    ) +
    legend_theme +
    tight_margin_mid

  # Panel C: pathway heatmap
  pathway_order <- intersect(pathway_order, rownames(pathway_z_mat))

  df_pw_long <- as.data.frame(pathway_z_mat[pathway_order, , drop = FALSE]) %>%
    tibble::rownames_to_column("pathway") %>%
    tidyr::pivot_longer(
      cols = -pathway,
      names_to = "bin",
      values_to = "z"
    ) %>%
    dplyr::mutate(
      bin = factor(as.integer(bin), levels = bin_levels),
      pathway_pretty = gsub("^(GOBP_|GO_BP_)", "", pathway),
      pathway_pretty = factor(
        pathway_pretty,
        levels = rev(gsub("^(GOBP_|GO_BP_)", "", pathway_order))
      )
    )

  pw_labs <- c(
    "CYTOPLASMIC_TRANSLATION" = "Cytoplasmic translation",
    "RIBOSOME_BIOGENESIS" = "Ribosome biogenesis",
    "RNA_SPLICING" = "RNA splicing",
    "POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_PROLIFERATION" = "αβ T-cell proliferation",
    "REGULATION_OF_T_CELL_ACTIVATION" = "T-cell activation",
    "LYMPHOCYTE_DIFFERENTIATION" = "Lymphocyte differentiation",
    "REGULATION_OF_LAMELLIPODIUM_ASSEMBLY" = "Lamellipodium assembly",
    "REGULATION_OF_IMMUNE_EFFECTOR_PROCESS" = "Immune effector process",
    "FC_RECEPTOR_MEDIATED_STIMULATORY_SIGNALING_PATHWAY" = "Fc receptor signaling",
    "T_CELL_MEDIATED_CYTOTOXICITY" = "T-cell cytotoxicity"
  )

  p_path <- ggplot2::ggplot(
    df_pw_long,
    ggplot2::aes(x = bin, y = pathway_pretty, fill = z)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = rd_bu_rev_w,
      limits = c(-clip_path_z, clip_path_z),
      values = scales::rescale(seq(-clip_path_z, clip_path_z, length.out = 100)),
      breaks = c(-clip_path_z, clip_path_z),
      oob = scales::squish,
      name = "Enrichment (z)"
    ) +
    ggplot2::scale_y_discrete(labels = pw_labs) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        barwidth = legend_barwidth,
        barheight = legend_barheight
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Functional association"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0, face = "bold", size = 10)
    ) +
    legend_theme +
    tight_margin_mid

  # Panel D: MW line
  p_mw <- plot_mw(
    mw_sum = mw_sum,
    bin_levels = bin_levels,
    title = "Median estimated MW (kDa)",
    hide_x = TRUE
  ) +
    tight_margin_bot

  (p_iso / p_markers / p_path / p_mw) +
    patchwork::plot_layout(
      heights = c(1.4, 1.1, 1.1, 1.0)
    ) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 12)
    )
}

# ============================================================
# 6) LOAD OBJECT
# ============================================================

message("Loading object: ", input_rds)
obj <- readRDS(input_rds)

if (!celltype_col %in% colnames(obj@meta.data)) {
  stop("Metadata column not found: ", celltype_col)
}

obj <- join_layer_if_possible(obj, assay = gene_assay, layer = gene_tpm_layer)
obj <- join_layer_if_possible(obj, assay = isoform_assay, layer = isoform_tpm_layer)

# ============================================================
# 7) COMPUTE CD45 MW AND ISOFORM EXPRESSION
# ============================================================

mw_res <- compute_cd45_mw(
  obj = obj,
  assay = isoform_assay,
  tpm_layer = isoform_tpm_layer,
  mw_col = mw_col,
  iso_levels = c("RO", "RB", "RA", "RAB", "RBC", "RABC")
)

obj <- mw_res$obj

# Reorder expression columns for downstream display.
cells_use <- rownames(obj@meta.data)[obj@meta.data[[celltype_col]] == celltype_val]

mw_named <- stats::setNames(
  obj@meta.data[[mw_col]],
  rownames(obj@meta.data)
)

cells_use <- intersect(
  cells_use,
  rownames(mw_res$iso_expr)[is.finite(mw_named[rownames(mw_res$iso_expr)])]
)

message("Selected ", celltype_val, " cells with finite inferred MW: ", length(cells_use))

if (length(cells_use) == 0) {
  stop("No cells available after cell type and MW filtering.")
}

iso_expr_use <- mw_res$iso_expr[cells_use, iso_order, drop = FALSE]

# ============================================================
# 8) BIN BY DOMINANT ISOFORM
# ============================================================

bin_res <- make_bins_by_dominant_isoform(
  iso_expr_cells_by_iso = iso_expr_use,
  mw_vec = mw_named,
  target_cells_per_bin = target_cells_per_bin,
  dom_order = iso_order,
  min_total_cd45_tpm = min_total_cd45_tpm,
  seed = seed_bins
)

bin_df <- bin_res$bin_df
iso_frac_per_cell <- bin_res$iso_frac_per_cell

message("Cells per bin:")
print(table(bin_df$bin))

message("Dominant isoform counts:")
print(table(bin_df$dominant_isoform))

message("NA MW in bin_df: ", sum(is.na(bin_df$MW)))

# ============================================================
# 9) BIN-LEVEL MW SUMMARY
# ============================================================

mw_sum <- bin_df %>%
  dplyr::group_by(bin) %>%
  dplyr::summarise(
    medMW = median(MW, na.rm = TRUE),
    meanMW = mean(MW, na.rm = TRUE),
    n_cells = dplyr::n(),
    dominant_isoform = as.character(dominant_isoform[1]),
    .groups = "drop"
  )

# ============================================================
# 10) BIN-LEVEL ISOFORM FRACTIONS
# ============================================================

iso_frac <- sapply(sort(unique(bin_df$bin)), function(b) {
  cells_b <- bin_df$cell[bin_df$bin == b]

  colMeans(
    iso_frac_per_cell[cells_b, iso_order, drop = FALSE],
    na.rm = TRUE
  )
})

iso_frac <- as.matrix(iso_frac)
rownames(iso_frac) <- iso_order
colnames(iso_frac) <- sort(unique(bin_df$bin))

iso_frac <- sweep(
  iso_frac,
  2,
  pmax(colSums(iso_frac), 1e-12),
  "/"
)

# ============================================================
# 11) AUCELL FUNCTIONAL ASSOCIATION
# ============================================================

cells_auc <- intersect(bin_df$cell, colnames(obj))
obj_sub <- subset(obj, cells = cells_auc)

gene_mat <- get_layer_safe(
  seu = obj_sub,
  assay = gene_assay,
  layer = gene_expr_layer,
  slot_fallback = gene_expr_layer
)

if (is.null(gene_mat)) {
  gene_mat <- get_layer_safe(
    seu = obj_sub,
    assay = "RNA",
    layer = gene_expr_layer,
    slot_fallback = gene_expr_layer
  )
}

if (is.null(gene_mat)) {
  stop("Could not read gene expression matrix for AUCell.")
}

# Retrieve GO:BP pathways from MSigDB through msigdbr.
go_bp <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "GO:BP"
)

all_pathways <- split(go_bp$gene_symbol, go_bp$gs_name)

missing_pw <- setdiff(pw_keep, names(all_pathways))
if (length(missing_pw) > 0) {
  warning("Requested pathways not found in MSigDB: ", paste(missing_pw, collapse = ", "))
}

pathways_list <- all_pathways[intersect(pw_keep, names(all_pathways))]

pathways_list <- lapply(pathways_list, function(x) {
  intersect(unique(x), rownames(gene_mat))
})

pathways_list <- pathways_list[lengths(pathways_list) >= 5]

if (length(pathways_list) == 0) {
  stop("No requested pathways have at least 5 genes present in the expression matrix.")
}

message("Pathways used for AUCell:")
print(names(pathways_list))

rankings <- tryCatch(
  AUCell::AUCell_buildRankings(
    gene_mat,
    nCores = 1,
    plotStats = FALSE
  ),
  error = function(e) {
    message("Retrying AUCell_buildRankings after converting to dense matrix.")
    AUCell::AUCell_buildRankings(
      as.matrix(gene_mat),
      nCores = 1,
      plotStats = FALSE
    )
  }
)

cellsAUC <- AUCell::AUCell_calcAUC(
  pathways_list,
  rankings,
  nCores = 1
)

auc_mat <- AUCell::getAUC(cellsAUC)

bin_df_auc <- bin_df %>%
  dplyr::filter(cell %in% colnames(auc_mat))

pw_mean <- summarise_auc_by_bin(
  auc_mat = auc_mat,
  bin_df = bin_df_auc
)

pathway_z_mat <- t(scale(t(pw_mean)))
pathway_z_mat[!is.finite(pathway_z_mat)] <- 0

# ============================================================
# 12) DISPLAY ORDER
# ============================================================

bin_levels_dom_then_mw <- mw_sum %>%
  dplyr::mutate(
    dominant_isoform = factor(dominant_isoform, levels = iso_order)
  ) %>%
  dplyr::arrange(dominant_isoform, dplyr::desc(medMW)) %>%
  dplyr::pull(bin)

# ============================================================
# 13) MAKE FIGURE
# ============================================================

plot_args <- list(
  iso_frac = iso_frac,
  mw_sum = mw_sum,
  pathway_z_mat = pathway_z_mat,
  obj = obj,
  bin_df = bin_df,
  marker_genes = marker_genes,
  clip_marker_z = clip_marker_z,
  clip_path_z = clip_path_z,
  iso_order = iso_order,
  iso_labels = iso_labels,
  pathway_order = names(pathways_list)
)

p_dom_order <- do.call(
  make_story_plot,
  c(
    list(bin_levels = bin_levels_dom_then_mw),
    plot_args
  )
)

print(p_dom_order)

# ============================================================
# 14) SAVE FIGURE AND SUPPORTING TABLES
# ============================================================

ggplot2::ggsave(
  filename = outfile_png,
  plot = p_dom_order,
  width = fig_width,
  height = fig_height,
  dpi = fig_dpi,
  bg = "white",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = outfile_pdf,
  plot = p_dom_order,
  width = fig_width,
  height = fig_height,
  device = cairo_pdf,
  bg = "white"
)

readr::write_csv(
  bin_df,
  file.path(output_dir, "Figure3_CD45_isoform_bins.csv")
)

readr::write_csv(
  mw_sum,
  file.path(output_dir, "Figure3_CD45_bin_MW_summary.csv")
)

readr::write_csv(
  as.data.frame(iso_frac) %>%
    tibble::rownames_to_column("isoform"),
  file.path(output_dir, "Figure3_CD45_bin_isoform_fractions.csv")
)

readr::write_csv(
  as.data.frame(pw_mean) %>%
    tibble::rownames_to_column("pathway"),
  file.path(output_dir, "Figure3_CD45_bin_pathway_AUC_mean.csv")
)

readr::write_csv(
  as.data.frame(pathway_z_mat) %>%
    tibble::rownames_to_column("pathway"),
  file.path(output_dir, "Figure3_CD45_bin_pathway_AUC_z.csv")
)

readr::write_csv(
  mw_res$iso_mass_tbl,
  file.path(output_dir, "Figure3_CD45_isoform_mass_table.csv")
)

readr::write_csv(
  mw_res$map_in_assay,
  file.path(output_dir, "Figure3_CD45_transcript_to_isoform_mapping_used.csv")
)

message("Saved: ", outfile_pdf)
message("Saved: ", outfile_png)
message("Done.")
