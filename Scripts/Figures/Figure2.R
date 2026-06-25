 #!/usr/bin/env Rscript

# Figure 2
# ============================================================
# FIGURE 2: Isoform diversity limits interpretation of
# CD marker protein-level measurements
#
# Panels:
# A = predicted localization of selected CD marker protein-coding transcripts
# B = localization-class switching across immune populations
# C = antibody PrEST coverage across PBMC-retained CD protein isoforms
#
# Saves:
# - Combined figure as PDF + PNG
# - Each panel separately as PDF + PNG
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
    library(ggrepel)
    library(purrr)
})


# ============================================================
# SETTINGS
# ============================================================

isoform_assay <- "Isoform"
counts_layer  <- "counts"
signal_layer  <- "tpm"

celltype_var <- "broad_celltype_final"

# Set to FALSE if you want no panel titles.
show_panel_titles <- TRUE

panel_title <- function(x) {
    if (isTRUE(show_panel_titles)) x else NULL
}

exclude_celltypes <- c("NK/ILC", "gdT", "Platelet", "MAIT", "pDC")

celltype_map <- c(
    "B_cell"   = "B cell",
    "CD4_T"    = "CD4 T",
    "CD8_T"    = "CD8 T",
    "DC"       = "DC",
    "Monocyte" = "Monocyte",
    "NK"       = "NK",
    "pDC"      = "pDC",
    "Treg"     = "Treg"
)

# Panel B raw order: pDC removed
ct_levels_switch_raw <- c(
    "CD4_T",
    "CD8_T",
    "Treg",
    "B_cell",
    "NK",
    "DC",
    "Monocyte"
)

# Panel B displayed order: pDC removed
ct_levels_switch_pretty <- c(
    "CD4 T",
    "CD8 T",
    "Treg",
    "B cell",
    "NK",
    "DC",
    "Monocyte"
)

ct_pretty_map <- c(
    "B_cell"   = "B cell",
    "CD4_T"    = "CD4 T",
    "CD8_T"    = "CD8 T",
    "DC"       = "DC",
    "Monocyte" = "Monocyte",
    "NK"       = "NK",
    "Treg"     = "Treg"
)

# Selected markers for Panel A
selected_cd_A <- c(
    "CD141", "CD161", "CD1c", "CD25", "CD79a", "CD66b", "CD2",
    "CD19", "CD14", "CD64", "CD10", "CD3e", "CD3d", "CD8a",
    "CD8b", "CD7", "CD303", "CD3g", "CD79b", "CD11c", "CD33",
    "CD11b", "CD20", "CD4", "CD56", "CD45", "CD123", "CD22",
    "CD247"
)

# ============================================================
# TEXT SIZES — MATCH FIGURE 1
# ============================================================

base_size <- 18

axis_text_size     <- 20
axis_title_size    <- 20
plot_title_size    <- 20
subtitle_size      <- 20
legend_text_size   <- 20
legend_title_size  <- 20
panel_tag_size     <- 22

# ggplot2 geom_text() sizes are not in points.
geom_text_size <- axis_text_size / ggplot2::.pt

# Panel-specific
panelA_axis_text_size   <- axis_text_size
panelA_legend_text_size <- legend_text_size

panelB_celltype_text_size <- axis_text_size
panelB_number_size        <- geom_text_size

# Scatter labels are kept slightly smaller than axis text to avoid crowding.
panelC_label_size <- 5.8

# ============================================================
# COLORS
# ============================================================

col_text <- "grey15"

col_membrane  <- "#74A9CF"
col_secreted  <- "#E08A66"
col_intracell <- "#E5E5B0"

col_full    <- "#2D6A4F"
col_partial <- "#B65C3A"

col_panelA <- "#4D4D4D"

# Grey scale matching Figure 1 isoform-switching triangle
col_grey_low  <- "#F2F2F2"
col_grey_high <- col_panelA

# ============================================================
# THEMES
# ============================================================

theme_fig <- ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
        text = ggplot2::element_text(color = col_text),
        axis.text = ggplot2::element_text(color = col_text, size = axis_text_size),
        axis.title = ggplot2::element_text(color = col_text, size = axis_title_size),
        plot.title = ggplot2::element_text(
            face = "plain",
            size = plot_title_size,
            color = col_text
        ),
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
        plot.title = ggplot2::element_text(
            face = "plain",
            size = plot_title_size,
            color = col_text
        ),
        plot.subtitle = ggplot2::element_text(size = subtitle_size, color = "grey35"),
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

standardize_loc_table <- function(loc_raw) {
    loc_raw %>%
        dplyr::mutate(
            tx_id = normalize_enst(transcript_id),
            compartment = as.character(compartment),
            prob = as.numeric(prob),
            loc_category = dplyr::case_when(
                compartment == "Cell membrane" ~ "Membrane",
                compartment == "Extracellular" ~ "Secreted",
                TRUE ~ "Intracellular"
            )
        ) %>%
        dplyr::filter(
            !is.na(tx_id),
            tx_id != "",
            !is.na(prob),
            !is.na(loc_category)
        ) %>%
        dplyr::group_by(tx_id, hgnc_symbol, loc_category) %>%
        dplyr::summarise(
            loc_prob = max(prob, na.rm = TRUE),
            top_compartment_in_class = compartment[which.max(prob)],
            .groups = "drop"
        ) %>%
        dplyr::group_by(tx_id) %>%
        dplyr::arrange(
            dplyr::desc(loc_prob),
            factor(loc_category, levels = c("Membrane", "Secreted", "Intracellular"))
        ) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select(
            tx_id,
            hgnc_symbol,
            loc_category,
            loc_prob,
            top_compartment_in_class
        )
}

# ============================================================
# LOAD INPUTS
# ============================================================

message("Loading transcript universe...")
tx_universe_pbmc <- readr::read_tsv(tx_universe_file, show_col_types = FALSE)

message("Loading CD marker metadata...")
cd_meta <- readr::read_csv(cd_marker_master_path, show_col_types = FALSE) %>%
    dplyr::transmute(
        cd_name = stringr::str_trim(`CD Marker`),
        gene_symbol = stringr::str_trim(`Gene Symbol`)
    ) %>%
    dplyr::distinct()

message("Loading localization table...")
loc_raw <- readRDS(loc_file)
loc_tbl <- standardize_loc_table(loc_raw)

message("Skipping Seurat object load; Figure 2 uses precomputed switching and antibody coverage files.")
# obj <- readRDS(seurat_file)

message("Loading antibody coverage files...")
ab_primary <- readRDS(ab_primary_file)
ab_gene_cov <- readRDS(ab_gene_cov_file)

# ============================================================
# COMMON TRANSCRIPT ANNOTATION
# ============================================================

tx_annot <- tx_universe_pbmc %>%
    dplyr::transmute(
        tx_id = normalize_enst(ensembl_transcript_id),
        gene_symbol = as.character(gene_symbol),
        transcript_biotype = as.character(transcript_biotype)
    ) %>%
    dplyr::filter(!is.na(tx_id), tx_id != "") %>%
    dplyr::distinct()

protein_coding_biotypes <- c(
    "protein_coding",
    "protein_coding_CDS_not_defined",
    "protein_coding_LoF"
)

tx_pc_loc <- tx_annot %>%
    dplyr::filter(transcript_biotype %in% protein_coding_biotypes) %>%
    dplyr::inner_join(loc_tbl, by = "tx_id") %>%
    dplyr::inner_join(cd_meta, by = "gene_symbol") %>%
    dplyr::distinct(tx_id, gene_symbol, cd_name, loc_category)

# ============================================================
# PANEL A — PREDICTED LOCALIZATION
# ============================================================

selected_meta_A <- cd_meta %>%
    dplyr::filter(cd_name %in% selected_cd_A) %>%
    dplyr::mutate(cd_name = factor(cd_name, levels = selected_cd_A))

plot_tbl_A <- tx_pc_loc %>%
    dplyr::inner_join(
        selected_meta_A %>% dplyr::select(gene_symbol, cd_name),
        by = c("gene_symbol", "cd_name")
    ) %>%
    dplyr::count(cd_name, gene_symbol, loc_category, name = "n_tx") %>%
    tidyr::complete(
        cd_name,
        gene_symbol,
        loc_category = c("Membrane", "Secreted", "Intracellular"),
        fill = list(n_tx = 0)
    ) %>%
    dplyr::group_by(cd_name, gene_symbol) %>%
    dplyr::mutate(total_tx = sum(n_tx)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(total_tx > 0)

order_A <- plot_tbl_A %>%
    dplyr::distinct(cd_name, total_tx) %>%
    dplyr::arrange(total_tx, cd_name) %>%
    dplyr::pull(cd_name) %>%
    as.character()

plot_tbl_A <- plot_tbl_A %>%
    dplyr::mutate(
        cd_name = factor(as.character(cd_name), levels = order_A),
        loc_category = factor(
            loc_category,
            levels = c("Intracellular", "Secreted", "Membrane")
        )
    )

pA <- ggplot2::ggplot(
    plot_tbl_A,
    ggplot2::aes(x = cd_name, y = n_tx, fill = loc_category)
) +
    ggplot2::geom_col(
        width = 0.84,
        color = "grey55",
        linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(
        values = c(
            "Membrane" = col_membrane,
            "Secreted" = col_secreted,
            "Intracellular" = col_intracell
        ),
        name = NULL
    ) +
    ggplot2::scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 1),
        breaks = c(1, 2, 3, 5, 10, 20, 30, 50, 100),
        limits = c(0, NA),
        expand = ggplot2::expansion(mult = c(0.02, 0.08))
    ) +
    ggplot2::labs(
        x = NULL,
        y = "Protein-coding transcripts",
        title = panel_title("Predicted localization of selected CD marker genes")
    ) +
    theme_fig +
    ggplot2::theme(
        plot.title = ggplot2::element_text(
            face = "plain",
            size = plot_title_size,
            color = col_text,
            margin = ggplot2::margin(t = 0, b = 10)
        ),
        axis.text.x = ggplot2::element_text(
            size = axis_text_size,
            angle = 45,
            hjust = 1,
            vjust = 1
        ),
        axis.text.y = ggplot2::element_text(size = axis_text_size),
        axis.title.x = ggplot2::element_text(size = axis_title_size),
        axis.title.y = ggplot2::element_text(size = axis_title_size),
        legend.text = ggplot2::element_text(size = legend_text_size),
        legend.title = ggplot2::element_text(size = legend_title_size),
        plot.margin = ggplot2::margin(18, 2, 44, 5.5)
    )

# ============================================================
# PANEL B — LOCALIZATION-CLASS SWITCHING TRIANGLE
# ============================================================

loc_switch_dir <- "/home/projects/dtu_00062/people/andson/SS3_PBMC/CD_markers_paper/Ensembl115/Results/Transcript_switching_analysis/Localization_class_switching_broad_celltype_final"

loc_switch_gene_file <- file.path(
    loc_switch_dir,
    "localization_class_all_top_switch_genes.rds"
)

loc_switch_pairs_file <- file.path(
    loc_switch_dir,
    "localization_class_all_switch_pairs.rds"
)

if (!file.exists(loc_switch_gene_file)) {
    stop("Missing localization-class switch results: ", loc_switch_gene_file)
}

loc_switch_genes <- readRDS(loc_switch_gene_file)

loc_switch_pairs <- if (file.exists(loc_switch_pairs_file)) {
    readRDS(loc_switch_pairs_file)
} else {
    NULL
}

pair_counts_obs <- loc_switch_genes %>%
    dplyr::filter(
        condition_1 %in% ct_levels_switch_raw,
        condition_2 %in% ct_levels_switch_raw
    ) %>%
    dplyr::mutate(
        i1 = match(condition_1, ct_levels_switch_raw),
        i2 = match(condition_2, ct_levels_switch_raw)
    ) %>%
    dplyr::filter(!is.na(i1), !is.na(i2), i1 != i2) %>%
    dplyr::mutate(
        ct_left_raw  = dplyr::if_else(i1 < i2, condition_1, condition_2),
        ct_right_raw = dplyr::if_else(i1 < i2, condition_2, condition_1)
    ) %>%
    dplyr::group_by(ct_left_raw, ct_right_raw) %>%
    dplyr::summarise(
        n_switching_genes = dplyr::n_distinct(gene_name),
        min_q = min(min_class_q, na.rm = TRUE),
        max_abs_dIF = max(max_abs_dIF, na.rm = TRUE),
        .groups = "drop"
    )

pair_grid <- purrr::map_dfr(
    combn(ct_levels_switch_raw, 2, simplify = FALSE),
    function(ct_pair) {
        tibble::tibble(
            ct_left_raw = ct_pair[1],
            ct_right_raw = ct_pair[2]
        )
    }
)

pair_counts <- pair_grid %>%
    dplyr::left_join(
        pair_counts_obs,
        by = c("ct_left_raw", "ct_right_raw")
    ) %>%
    dplyr::mutate(
        n_switching_genes = dplyr::coalesce(n_switching_genes, 0L),
        ct_left = dplyr::recode(ct_left_raw, !!!ct_pretty_map),
        ct_right = dplyr::recode(ct_right_raw, !!!ct_pretty_map),
        x = factor(ct_left, levels = ct_levels_switch_pretty),
        y = factor(ct_right, levels = rev(ct_levels_switch_pretty)),
        label = dplyr::if_else(
            n_switching_genes > 0,
            as.character(n_switching_genes),
            ""
        ),
        fill_value = dplyr::if_else(
            n_switching_genes > 0,
            as.numeric(n_switching_genes),
            NA_real_
        )
    )

max_count_B <- max(pair_counts$n_switching_genes, na.rm = TRUE)

pair_counts <- pair_counts %>%
    dplyr::mutate(
        label_col = dplyr::if_else(
            max_count_B > 0 & n_switching_genes >= 0.60 * max_count_B,
            "white",
            "grey15"
        )
    )

pB <- ggplot2::ggplot(
    pair_counts,
    ggplot2::aes(x = x, y = y, fill = fill_value)
) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
        ggplot2::aes(label = label, color = label_col),
        size = panelB_number_size
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient(
        low = col_grey_low,
        high = col_grey_high,
        na.value = "white",
        name = "# genes with\nlocalization-class\nswitch",
        guide = ggplot2::guide_colorbar(
            barheight = grid::unit(52, "mm"),
            barwidth = grid::unit(8, "mm"),
            title.position = "top",
            title.hjust = 0.5,
            label.position = "right"
        )
    ) +
    ggplot2::labs(
        x = NULL,
        y = NULL,
        title = panel_title("Localization-class switching across immune populations")
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
            angle = 60,
            hjust = 1,
            vjust = 1,
            color = col_text,
            size = panelB_celltype_text_size
        ),
        axis.text.y = ggplot2::element_text(
            color = col_text,
            size = panelB_celltype_text_size
        ),
        plot.title = ggplot2::element_text(
            face = "plain",
            size = plot_title_size,
            color = col_text,
            margin = ggplot2::margin(t = 0, b = 8)
        ),
        plot.subtitle = ggplot2::element_text(
            size = subtitle_size,
            color = "grey35"
        ),
        legend.title = ggplot2::element_text(size = legend_title_size),
        legend.text = ggplot2::element_text(size = legend_text_size),
        plot.margin = ggplot2::margin(18, 24, 5.5, 5.5)
    ) +
    ggplot2::coord_equal()

# ============================================================
# PANEL C — ANTIBODY PrEST COVERAGE SCATTER
# ============================================================

plot_df_C <- ab_primary %>%
    dplyr::filter(n_iso_hit > 0) %>%
    dplyr::left_join(cd_meta, by = "gene_symbol") %>%
    dplyr::mutate(
        marker_label = dplyr::if_else(
            !is.na(cd_name) & cd_name != "",
            cd_name,
            gene_symbol
        ),
        coverage_class = dplyr::if_else(
            frac_hit_iso == 1,
            "Full coverage",
            "Partial coverage"
        ),
        coverage_class = factor(
            coverage_class,
            levels = c("Full coverage", "Partial coverage")
        )
    ) %>%
    dplyr::filter(
        !is.na(ab_id),
        !is.na(gene_symbol),
        !is.na(n_iso_total),
        !is.na(frac_hit_iso),
        n_iso_total > 0
    )

# Only label the rightmost/high-complexity markers.
rightmost_labels <- plot_df_C %>%
    dplyr::filter(
        n_iso_total >= 20,
        frac_hit_iso < 0.95
    ) %>%
    dplyr::mutate(label_score = n_iso_total * (1 - frac_hit_iso)) %>%
    dplyr::arrange(desc(n_iso_total), desc(label_score)) %>%
    dplyr::slice_head(n = 10)

# CD100 and CD88 are intentionally not forced/labeled.
force_cds_right <- c(
    "CD44",
    "CD99",
    "CD224",
    "CD36",
    "CD82",
    "CD123",
    "CD247",
    "CD55"
)

label_df_C <- dplyr::bind_rows(
    rightmost_labels,
    plot_df_C %>% dplyr::filter(marker_label %in% force_cds_right)
) %>%
    dplyr::distinct(marker_label, .keep_all = TRUE) %>%
    dplyr::filter(!marker_label %in% c("CD100", "CD88"))

pC <- ggplot2::ggplot(
    plot_df_C,
    ggplot2::aes(x = n_iso_total, y = frac_hit_iso)
) +
    ggplot2::geom_jitter(
        ggplot2::aes(color = coverage_class),
        width = 0.08,
        height = 0.015,
        size = 2.2,
        alpha = 0.75
    ) +
    ggrepel::geom_text_repel(
        data = label_df_C,
        ggplot2::aes(label = marker_label),
        size = panelC_label_size,
        color = col_text,
        segment.color = NA,
        box.padding = 0.35,
        point.padding = 0.18,
        max.overlaps = Inf,
        seed = 1
    ) +
    ggplot2::scale_color_manual(
        values = c(
            "Full coverage" = col_full,
            "Partial coverage" = col_partial
        ),
        name = NULL
    ) +
    ggplot2::scale_x_continuous(
        trans = "log10",
        breaks = c(1, 2, 3, 5, 10, 20, 50, 100)
    ) +
    ggplot2::scale_y_continuous(
        limits = c(-0.04, 1.05),
        breaks = seq(0, 1, 0.25),
        labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
        x = "PBMC protein isoforms for antibody target",
        y = "Target isoforms",
        title = panel_title("PrEST coverage across CD protein isoforms")
    ) +
    theme_fig +
    ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(size = axis_text_size),
        axis.text.y = ggplot2::element_text(size = axis_text_size),
        axis.title.x = ggplot2::element_text(size = axis_title_size),
        axis.title.y = ggplot2::element_text(
            size = axis_title_size,
            margin = ggplot2::margin(r = 14)
        ),
        plot.title = ggplot2::element_text(
            face = "plain",
            size = plot_title_size,
            color = col_text,
            margin = ggplot2::margin(t = 0, b = 8)
        ),
        legend.text = ggplot2::element_text(size = legend_text_size),
        legend.title = ggplot2::element_text(size = legend_title_size),
        plot.margin = ggplot2::margin(18, 5.5, 5.5, 24)
    )

# ============================================================
# SAVE SEPARATE PANELS
# ============================================================

save_panel(
    pA,
    "Figure2A_predicted_localization_selected_CD_markers",
    width = 18.0,
    height = 4.8,
    dpi = 600
)

save_panel(
    pB,
    "Figure2B_localization_class_switching",
    width = 8.4,
    height = 6.6,
    dpi = 600
)

save_panel(
    pC,
    "Figure2C_antibody_PrEST_coverage",
    width = 9.4,
    height = 6.6,
    dpi = 600
)

# ============================================================
# COMBINED FIGURE
# ============================================================

bottom_row <- (pB | pC) +
    patchwork::plot_layout(widths = c(1.05, 1.35))

pA_for_combined <- patchwork::free(pA, side = "r")

fig_all <- pA_for_combined / bottom_row +
    patchwork::plot_layout(
        heights = c(1.12, 1.88)
    ) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(
        plot.tag = ggplot2::element_text(
            face = "bold",
            size = panel_tag_size,
            color = col_text,
            vjust = 0.5
        ),
        plot.tag.position = c(0, 0.96),
        plot.margin = ggplot2::margin(5.5, 10, 5.5, 5.5)
    )

print(fig_all)

ggplot2::ggsave(
    outfile_combined_pdf,
    fig_all,
    width = 18.0,
    height = 16.8,
    device = cairo_pdf,
    bg = "white"
)

ggplot2::ggsave(
    outfile_combined_png,
    fig_all,
    width = 18.0,
    height = 16.8,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
)

message("Saved combined PDF: ", outfile_combined_pdf)
message("Saved combined PNG: ", outfile_combined_png)
