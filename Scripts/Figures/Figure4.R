# ============================================================
# Figure4.R
# IL7R/CD127 + CD74 isoform landscape, protein structure,
# PrEST mapping, and functional programs
#
#
# Required input files:
#   data/processed/figure4_inputs/tx_summary.rds
#   data/processed/figure4_inputs/function_main_bundle.rds
#   data/processed/figure4_inputs/tm_seg2.rds
#   data/processed/figure4_inputs/ipr_df_pfam.rds
#   data/processed/figure4_inputs/ab_tx.rds
#   data/processed/figure4_inputs/cd_markers_localization.rds
#
# Output:
#   results/figures/figure4/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(cowplot)
  library(purrr)
  library(msa)
  library(Biostrings)
  library(ggnewscale)
  library(grid)
})

# ============================================================
# PATHS
# ============================================================

input_dir <- file.path("data", "processed", "figure4_inputs")
fig_out_dir <- file.path("results", "figures", "figure4")
structure_dir <- file.path(fig_out_dir, "prest_structure_examples")

dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(structure_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- list(
  tx_summary = file.path(input_dir, "tx_summary.rds"),
  function_main_bundle = file.path(input_dir, "function_main_bundle.rds"),
  tm_seg2 = file.path(input_dir, "tm_seg2.rds"),
  ipr_df_pfam = file.path(input_dir, "ipr_df_pfam.rds"),
  ab_tx = file.path(input_dir, "ab_tx.rds"),
  cd_markers_localization = file.path(input_dir, "cd_markers_localization.rds")
)

missing_input_files <- input_files[!vapply(input_files, file.exists, logical(1))]

if (length(missing_input_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(unlist(missing_input_files), collapse = "\n"),
    "\nRun 01_make_figure4_input_objects.R first."
  )
}

# ============================================================
# SETTINGS
# ============================================================

cell_order <- c(
  "CD4 T",
  "CD8 T",
  "Treg",
  "B cell",
  "NK",
  "DC",
  "Monocyte"
)

text_col <- "grey30"
col_other <- "#6B6B6B"

nes_low  <- "#3B6FB6"
nes_mid  <- "white"
nes_high <- "#B22222"

TEXT_SIZE <- 14
GEOM_TEXT_SIZE <- TEXT_SIZE / ggplot2::.pt
TITLE_SIZE <- 17
MAIN_TITLE_SIZE <- 19
TAG_SIZE <- 14

base_size <- TEXT_SIZE

top_cols <- c(
  "#2D6A4F",
  "#4C78A8",
  "#F58518",
  "#72B7B2",
  "#B65C3A",
  "#8E6C8A",
  "#54A24B",
  "#E45756",
  "#A0CBE8",
  "#FFBE7D"
)

loc_cols <- c(
  "Membrane"      = "#74A9CF",
  "Secreted"      = "#E08A66",
  "Intracellular" = "#E5E5B0",
  "No prediction" = "grey85"
)

# ============================================================
# BASIC HELPERS
# ============================================================

strip_tx_version <- function(x) {
  sub("\\..*$", "", as.character(x))
}

extract_enst <- function(x) {
  out <- stringr::str_extract(as.character(x), "ENST\\d+")
  ifelse(is.na(out), as.character(x), out)
}

clean_aa <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- gsub("\\*", "", x)
  x <- gsub("-", "", x)
  x <- gsub("[^A-Z]", "", x)
  dplyr::na_if(x, "")
}

prettify_pathway <- function(x) {
  x %>%
    stringr::str_replace("^GOBP_", "") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("\\bil 7\\b", "IL-7") %>%
    stringr::str_replace_all("\\bmhc\\b", "MHC") %>%
    stringr::str_replace_all("\\bii\\b", "II") %>%
    stringr::str_replace_all("\\bt cell\\b", "T cell") %>%
    stringr::str_replace_all("\\bb cell\\b", "B cell") %>%
    stringr::str_to_sentence()
}

save_plot_pair <- function(plot, filename_base, width, height, dpi = 600) {
  pdf_file <- file.path(fig_out_dir, paste0(filename_base, ".pdf"))
  png_file <- file.path(fig_out_dir, paste0(filename_base, ".png"))
  
  ggplot2::ggsave(
    pdf_file,
    plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
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
# LOAD INPUT OBJECTS
# ============================================================

message("Loading Figure 4 input objects...")

tx_summary <- readRDS(input_files$tx_summary)
function_main_bundle <- readRDS(input_files$function_main_bundle)
tm_seg2 <- readRDS(input_files$tm_seg2)
ipr_df_pfam <- readRDS(input_files$ipr_df_pfam)
ab_tx <- readRDS(input_files$ab_tx)
cd_markers_localization <- readRDS(input_files$cd_markers_localization)

# ============================================================
# VALIDATE INPUT OBJECTS
# ============================================================

required_tx_summary_cols <- c(
  "transcript_id",
  "gene_symbol",
  "cell_type_display"
)

missing_tx_summary_cols <- setdiff(required_tx_summary_cols, colnames(tx_summary))

if (length(missing_tx_summary_cols) > 0) {
  stop(
    "tx_summary is missing required column(s): ",
    paste(missing_tx_summary_cols, collapse = ", ")
  )
}

if (!"total_signal" %in% colnames(tx_summary)) {
  signal_candidates <- c(
    "total_signal",
    "tx_signal",
    "transcript_signal",
    "signal",
    "sum_signal",
    "sum_tpm",
    "total_tpm",
    "mean_tpm",
    "tpm",
    "counts",
    "count",
    "n_counts",
    "abundance",
    "isoform_fraction"
  )
  
  signal_col <- intersect(signal_candidates, colnames(tx_summary))[1]
  
  if (is.na(signal_col)) {
    stop(
      "tx_summary needs a total_signal column or one usable signal column. Available columns are:\n",
      paste(colnames(tx_summary), collapse = ", ")
    )
  }
  
  message("Creating tx_summary$total_signal from column: ", signal_col)
  
  tx_summary <- tx_summary %>%
    dplyr::mutate(total_signal = as.numeric(.data[[signal_col]]))
}

required_function_bundle_cols <- c(
  "gene",
  "population",
  "isoform",
  "pathway",
  "NES",
  "padj",
  "is_main_pathway"
)

missing_function_bundle_cols <- setdiff(
  required_function_bundle_cols,
  colnames(function_main_bundle)
)

if (length(missing_function_bundle_cols) > 0) {
  stop(
    "function_main_bundle is missing required column(s): ",
    paste(missing_function_bundle_cols, collapse = ", ")
  )
}

required_tm_seg2_cols <- c(
  "transcript_id",
  "hgnc_symbol",
  "state",
  "start",
  "end"
)

missing_tm_seg2_cols <- setdiff(required_tm_seg2_cols, colnames(tm_seg2))

if (length(missing_tm_seg2_cols) > 0) {
  stop(
    "tm_seg2 is missing required column(s): ",
    paste(missing_tm_seg2_cols, collapse = ", ")
  )
}

required_pfam_cols <- c(
  "transcript_id",
  "domain_label",
  "start",
  "end"
)

missing_pfam_cols <- setdiff(required_pfam_cols, colnames(ipr_df_pfam))

if (length(missing_pfam_cols) > 0) {
  warning(
    "ipr_df_pfam is missing Pfam column(s): ",
    paste(missing_pfam_cols, collapse = ", "),
    ". Pfam domains will not be drawn."
  )
}

required_ab_tx_cols <- c(
  "transcript_id",
  "gene_symbol",
  "ab_id",
  "protein_seq",
  "prest_seq",
  "present"
)

missing_ab_tx_cols <- setdiff(required_ab_tx_cols, colnames(ab_tx))

if (length(missing_ab_tx_cols) > 0) {
  stop(
    "ab_tx is missing required column(s): ",
    paste(missing_ab_tx_cols, collapse = ", ")
  )
}

required_loc_cols <- c(
  "transcript_id",
  "hgnc_symbol",
  "compartment",
  "prob"
)

missing_loc_cols <- setdiff(required_loc_cols, colnames(cd_markers_localization))

if (length(missing_loc_cols) > 0) {
  stop(
    "cd_markers_localization is missing required column(s): ",
    paste(missing_loc_cols, collapse = ", ")
  )
}

# ============================================================
# CLEAN INPUT OBJECTS
# ============================================================

tx_summary <- tx_summary %>%
  dplyr::mutate(
    transcript_id = strip_tx_version(transcript_id),
    gene_symbol = as.character(gene_symbol),
    cell_type_display = as.character(cell_type_display),
    total_signal = as.numeric(total_signal)
  ) %>%
  dplyr::filter(
    gene_symbol %in% c("IL7R", "CD74"),
    cell_type_display %in% cell_order,
    !is.na(transcript_id),
    transcript_id != "",
    !is.na(total_signal),
    total_signal > 0
  ) %>%
  dplyr::mutate(
    cell_type_display = factor(cell_type_display, levels = cell_order)
  )

function_main_bundle <- function_main_bundle %>%
  dplyr::mutate(
    gene = as.character(gene),
    population = as.character(population),
    isoform = strip_tx_version(isoform)
  )

tm_seg2 <- tm_seg2 %>%
  dplyr::mutate(
    transcript_id = strip_tx_version(transcript_id),
    hgnc_symbol = as.character(hgnc_symbol),
    state = as.character(state),
    start = as.integer(start),
    end = as.integer(end)
  )

ipr_df_pfam <- ipr_df_pfam %>%
  dplyr::mutate(
    transcript_id = strip_tx_version(transcript_id)
  )

ab_tx <- ab_tx %>%
  dplyr::mutate(
    transcript_id = strip_tx_version(transcript_id),
    gene_symbol = as.character(gene_symbol),
    ab_id = as.character(ab_id)
  )

if (!is.logical(ab_tx$present)) {
  ab_tx <- ab_tx %>%
    dplyr::mutate(
      present = as.character(present) %in% c("TRUE", "True", "T", "1")
    )
}

cd_markers_localization <- cd_markers_localization %>%
  dplyr::mutate(
    transcript_id = strip_tx_version(transcript_id),
    hgnc_symbol = as.character(hgnc_symbol),
    compartment = as.character(compartment),
    prob = as.numeric(prob)
  )

message("Panel A cell types retained:")
print(unique(tx_summary$cell_type_display))

# ============================================================
# LABEL MAPS
# ============================================================

tx_label_map_il7r_all <- c(
  "ENST00000303115" = "ENST00000303115",
  "ENST00000506850" = "ENST00000506850",
  "ENST00000511982" = "ENST00000511982",
  "ENST00000515665" = "ENST00000515665",
  "ENST00000514217" = "ENST00000514217"
)

tx_label_map_cd74_all <- c(
  "ENST00000881215" = "ENST00000881215",
  "ENST00000881218" = "ENST00000881218",
  "ENST00000353334" = "ENST00000353334",
  "ENST00000009530" = "ENST00000009530",
  "ENST00000377795" = "ENST00000377795"
)

make_tx_label_map_for_selected <- function(selected_tx, full_label_map) {
  selected_tx <- strip_tx_version(selected_tx)
  out <- setNames(selected_tx, selected_tx)
  overlapping <- intersect(selected_tx, names(full_label_map))
  out[overlapping] <- full_label_map[overlapping]
  out
}

# ============================================================
# TOP TRANSCRIPT HELPERS
# ============================================================

get_top_tx_for_gene_celltype <- function(tx_summary,
                                         target_gene,
                                         focal_celltype,
                                         top_n = 3) {
  tx_summary %>%
    dplyr::filter(
      gene_symbol == target_gene,
      cell_type_display == focal_celltype
    ) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(
      signal = sum(total_signal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(signal)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(transcript_id)
}

get_top_tx_per_celltype <- function(tx_summary,
                                    target_gene,
                                    top_n = 3) {
  tx_summary %>%
    dplyr::filter(
      gene_symbol == target_gene,
      cell_type_display %in% cell_order
    ) %>%
    dplyr::group_by(cell_type_display, transcript_id) %>%
    dplyr::summarise(
      signal = sum(total_signal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(cell_type_display) %>%
    dplyr::slice_max(
      order_by = signal,
      n = top_n,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(cell_type_display, transcript_id)
}

make_color_vector_for_landscape <- function(top_tx_tbl,
                                            full_label_map,
                                            other_col = "#6B6B6B") {
  all_top_tx <- unique(top_tx_tbl$transcript_id)
  
  label_map <- make_tx_label_map_for_selected(
    selected_tx = all_top_tx,
    full_label_map = full_label_map
  )
  
  tx_labels <- unname(label_map[all_top_tx])
  
  if (length(tx_labels) <= length(top_cols)) {
    cols <- top_cols[seq_along(tx_labels)]
  } else {
    cols <- scales::hue_pal()(length(tx_labels))
  }
  
  out <- setNames(cols, tx_labels)
  out <- c(out, "Other retained isoforms" = other_col)
  
  list(
    label_map = label_map,
    fill_cols = out
  )
}

# ============================================================
# PANEL A: TRANSCRIPT LANDSCAPE
# ============================================================

make_transcript_landscape_top3_per_celltype <- function(tx_summary,
                                                        target_gene,
                                                        gene_title,
                                                        full_label_map,
                                                        top_n = 3,
                                                        cell_order = NULL,
                                                        other_col = "#6B6B6B") {
  
  top_tx_tbl <- get_top_tx_per_celltype(
    tx_summary = tx_summary,
    target_gene = target_gene,
    top_n = top_n
  )
  
  color_info <- make_color_vector_for_landscape(
    top_tx_tbl = top_tx_tbl,
    full_label_map = full_label_map,
    other_col = other_col
  )
  
  label_map <- color_info$label_map
  fill_cols <- color_info$fill_cols
  
  plot_df <- tx_summary %>%
    dplyr::filter(
      gene_symbol == target_gene,
      !is.na(cell_type_display),
      cell_type_display %in% cell_order,
      !is.na(transcript_id),
      transcript_id != "",
      !is.na(total_signal),
      total_signal > 0
    ) %>%
    dplyr::group_by(cell_type_display, transcript_id) %>%
    dplyr::summarise(
      signal = sum(total_signal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      top_tx_tbl %>% dplyr::mutate(is_top = TRUE),
      by = c("cell_type_display", "transcript_id")
    ) %>%
    dplyr::mutate(
      is_top = dplyr::coalesce(is_top, FALSE),
      tx_label = dplyr::recode(
        transcript_id,
        !!!label_map,
        .default = transcript_id
      ),
      tx_group = dplyr::if_else(
        is_top,
        tx_label,
        "Other retained isoforms"
      )
    ) %>%
    dplyr::group_by(cell_type_display, tx_group) %>%
    dplyr::summarise(
      signal = sum(signal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(cell_type_display) %>%
    dplyr::mutate(
      plotted_total_signal = sum(signal, na.rm = TRUE),
      frac = dplyr::if_else(
        plotted_total_signal > 0,
        signal / plotted_total_signal,
        NA_real_
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(frac))
  
  if (!is.null(cell_order)) {
    plot_df <- plot_df %>%
      dplyr::mutate(
        cell_type_display = factor(cell_type_display, levels = cell_order)
      ) %>%
      dplyr::filter(!is.na(cell_type_display))
  }
  
  tx_order <- plot_df %>%
    dplyr::filter(tx_group != "Other retained isoforms") %>%
    dplyr::group_by(tx_group) %>%
    dplyr::summarise(total_frac = sum(frac, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total_frac)) %>%
    dplyr::pull(tx_group)
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      tx_group = factor(
        tx_group,
        levels = c(tx_order, "Other retained isoforms")
      )
    )
  
  message("Panel A fraction check for ", target_gene, ":")
  print(
    plot_df %>%
      dplyr::group_by(cell_type_display) %>%
      dplyr::summarise(
        plotted_sum = sum(frac, na.rm = TRUE),
        .groups = "drop"
      )
  )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = cell_type_display, y = frac, fill = tx_group)
  ) +
    ggplot2::geom_col(
      width = 0.78,
      color = "white",
      linewidth = 0.15
    ) +
    ggplot2::scale_fill_manual(
      values = fill_cols,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = gene_title,
      x = NULL,
      y = "Transcript fraction\nof gene signal",
      fill = "Transcript isoform"
    ) +
    ggplot2::theme_classic(base_size = TEXT_SIZE) +
    ggplot2::theme(
      text = ggplot2::element_text(color = text_col, size = TEXT_SIZE),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = TITLE_SIZE,
        hjust = 0.5,
        color = "grey20"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 35,
        hjust = 1,
        color = text_col,
        size = TEXT_SIZE
      ),
      axis.text.y = ggplot2::element_text(color = text_col, size = TEXT_SIZE),
      axis.title.y = ggplot2::element_text(color = text_col, size = TEXT_SIZE),
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(size = TEXT_SIZE),
      legend.text = ggplot2::element_text(
        size = TEXT_SIZE,
        margin = ggplot2::margin(b = 4)
      ),
      legend.key.size = grid::unit(0.48, "cm"),
      legend.key.height = grid::unit(0.56, "cm"),
      legend.spacing.y = grid::unit(0.16, "cm"),
      plot.margin = ggplot2::margin(4, 8, 4, 4)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        ncol = 1,
        byrow = TRUE,
        title.position = "top",
        keyheight = grid::unit(0.56, "cm")
      )
    )
}

# ============================================================
# MSA / STRUCTURE HELPERS
# ============================================================

msa_to_aastringset <- function(msa_obj, expected_names = NULL) {
  aln <- msa::msaConvert(msa_obj, type = "seqinr::alignment")
  aa <- Biostrings::AAStringSet(aln$seq)
  names(aa) <- strip_tx_version(extract_enst(aln$nam))
  
  if (!is.null(expected_names)) {
    expected_names <- strip_tx_version(expected_names)
    
    if (!all(expected_names %in% names(aa)) && length(aa) == length(expected_names)) {
      names(aa) <- expected_names
    }
    
    aa <- aa[expected_names]
  }
  
  aa
}

build_orig_to_aln_map <- function(aligned_seq_char) {
  chars <- strsplit(aligned_seq_char, "")[[1]]
  map <- integer(0)
  orig_pos <- 0L
  
  for (aln_pos in seq_along(chars)) {
    if (chars[aln_pos] != "-") {
      orig_pos <- orig_pos + 1L
      map[orig_pos] <- aln_pos
    }
  }
  
  map
}

build_topology_codes <- function(seg_df, pep_len) {
  codes <- rep("U", pep_len)
  
  priority <- c(outside = 1, inside = 1, signal = 2, TMhelix = 3)
  state_to_code <- c(signal = "S", TMhelix = "M", inside = "I", outside = "O")
  
  if (nrow(seg_df) == 0) {
    return(codes)
  }
  
  seg_df2 <- seg_df %>%
    dplyr::mutate(state = as.character(state)) %>%
    dplyr::filter(!is.na(start), !is.na(end), !is.na(state)) %>%
    dplyr::mutate(
      start = pmax(1L, as.integer(start)),
      end = pmin(as.integer(pep_len), as.integer(end)),
      pr = as.integer(priority[state])
    ) %>%
    dplyr::filter(!is.na(pr), start <= end) %>%
    dplyr::arrange(pr)
  
  if (nrow(seg_df2) == 0) {
    return(codes)
  }
  
  for (i in seq_len(nrow(seg_df2))) {
    st <- seg_df2$state[i]
    if (!st %in% names(state_to_code)) next
    
    s <- seg_df2$start[i]
    e <- seg_df2$end[i]
    
    if (is.na(s) || is.na(e) || s > e) next
    
    codes[s:e] <- state_to_code[[st]]
  }
  
  codes
}

map_codes_to_alignment <- function(aligned_seq_char, codes) {
  chars <- strsplit(aligned_seq_char, "")[[1]]
  out <- character(length(chars))
  orig_pos <- 1L
  
  for (aln_pos in seq_along(chars)) {
    if (chars[aln_pos] == "-") {
      out[aln_pos] <- "-"
    } else {
      out[aln_pos] <- if (orig_pos <= length(codes)) codes[orig_pos] else "U"
      orig_pos <- orig_pos + 1L
    }
  }
  
  out
}

code_to_type <- function(code) {
  dplyr::case_when(
    code == "S" ~ "Signal peptide",
    code == "O" ~ "Extracellular",
    code == "M" ~ "Transmembrane",
    code == "I" ~ "Intracellular",
    code == "-" ~ "Alignment gap",
    TRUE ~ "Unknown"
  )
}

find_prest_interval <- function(prest_seq, protein_seq, max_mismatch_frac = 0.05) {
  prest <- clean_aa(prest_seq)
  prot <- clean_aa(protein_seq)
  
  if (is.na(prest) || is.na(prot) || prest == "" || prot == "") {
    return(tibble::tibble(aa_start = NA_integer_, aa_end = NA_integer_))
  }
  
  hit <- regexpr(prest, prot, fixed = TRUE)[1]
  
  if (!is.na(hit) && hit > 0) {
    return(tibble::tibble(
      aa_start = as.integer(hit),
      aa_end = as.integer(hit + nchar(prest) - 1)
    ))
  }
  
  max_mismatch <- floor(max_mismatch_frac * nchar(prest))
  
  mp <- tryCatch(
    Biostrings::matchPattern(
      pattern = Biostrings::AAString(prest),
      subject = Biostrings::AAString(prot),
      max.mismatch = max_mismatch,
      fixed = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(mp) && length(mp) > 0) {
    return(tibble::tibble(
      aa_start = as.integer(Biostrings::start(mp)[1]),
      aa_end = as.integer(Biostrings::end(mp)[1])
    ))
  }
  
  tibble::tibble(aa_start = NA_integer_, aa_end = NA_integer_)
}

project_domains_to_alignment <- function(dom_raw,
                                         aligned,
                                         y_base,
                                         ymin_off = -0.44,
                                         ymax_off = -0.25) {
  if (nrow(dom_raw) == 0) {
    return(tibble::tibble())
  }
  
  purrr::pmap_dfr(dom_raw, function(transcript_id, domain_label, aa_start, aa_end, ...) {
    tid <- as.character(transcript_id)
    
    if (!tid %in% names(aligned)) return(tibble::tibble())
    
    aln_seq <- as.character(aligned[[tid]])
    map <- build_orig_to_aln_map(aln_seq)
    
    if (length(map) == 0) return(tibble::tibble())
    
    ungapped_len <- length(map)
    
    s <- max(1L, as.integer(aa_start))
    e <- min(as.integer(aa_end), ungapped_len)
    
    if (is.na(s) || is.na(e) || s > e) return(tibble::tibble())
    
    tibble::tibble(
      transcript_id = tid,
      domain_label = as.character(domain_label),
      aln_start = as.integer(map[s]),
      aln_end = as.integer(map[e]),
      ymin = y_base[[tid]] + ymin_off,
      ymax = y_base[[tid]] + ymax_off
    )
  }) %>%
    dplyr::distinct()
}

project_prest_to_alignment <- function(prest_raw,
                                       aligned,
                                       y_base,
                                       ymin_off = -0.24,
                                       ymax_off = 0.24) {
  if (nrow(prest_raw) == 0) {
    return(tibble::tibble())
  }
  
  purrr::pmap_dfr(prest_raw, function(transcript_id, aa_start, aa_end, present, ...) {
    tid <- as.character(transcript_id)
    
    if (!isTRUE(present)) return(tibble::tibble())
    if (!tid %in% names(aligned)) return(tibble::tibble())
    if (is.na(aa_start) || is.na(aa_end) || aa_start > aa_end) return(tibble::tibble())
    
    aln_seq <- as.character(aligned[[tid]])
    map <- build_orig_to_aln_map(aln_seq)
    
    if (length(map) == 0) return(tibble::tibble())
    
    ungapped_len <- length(map)
    
    s <- max(1L, as.integer(aa_start))
    e <- min(as.integer(aa_end), ungapped_len)
    
    if (is.na(s) || is.na(e) || s > e) return(tibble::tibble())
    
    aln_positions <- as.integer(map[s:e])
    
    if (length(aln_positions) == 0) return(tibble::tibble())
    
    grp <- cumsum(c(1, diff(aln_positions) != 1))
    
    tibble::tibble(
      transcript_id = tid,
      aln_pos = aln_positions,
      grp = grp
    ) %>%
      dplyr::group_by(transcript_id, grp) %>%
      dplyr::summarise(
        aln_start = min(aln_pos),
        aln_end = max(aln_pos),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        ymin = y_base[[tid]] + ymin_off,
        ymax = y_base[[tid]] + ymax_off,
        prest_status = "HPA PrEST"
      )
  }) %>%
    dplyr::distinct()
}

# ============================================================
# PANEL B: STRUCTURE / PrEST PLOT FUNCTION
# ============================================================

plot_hpa_prest_structure <- function(target_gene,
                                     target_cd,
                                     target_ab,
                                     ab_tx,
                                     tm_seg2,
                                     pfam_data,
                                     output_file_pdf,
                                     selected_tx_keep,
                                     max_mismatch_frac = 0.05) {
  
  selected_tx_keep <- strip_tx_version(selected_tx_keep)
  
  ab_sel <- ab_tx %>%
    dplyr::mutate(
      transcript_id = strip_tx_version(transcript_id),
      gene_symbol = as.character(gene_symbol),
      ab_id = as.character(ab_id),
      protein_seq = clean_aa(protein_seq),
      prest_seq = clean_aa(prest_seq)
    ) %>%
    dplyr::filter(
      gene_symbol == target_gene,
      ab_id == target_ab,
      transcript_id %in% selected_tx_keep,
      !is.na(transcript_id),
      transcript_id != "",
      !is.na(protein_seq),
      protein_seq != ""
    ) %>%
    dplyr::distinct(
      ab_id,
      gene_symbol,
      prest_seq,
      transcript_id,
      protein_isoform_id,
      protein_seq,
      present,
      .keep_all = TRUE
    ) %>%
    dplyr::mutate(
      tx_order_factor = factor(transcript_id, levels = selected_tx_keep),
      pep_len = nchar(protein_seq)
    ) %>%
    dplyr::arrange(tx_order_factor)
  
  if (nrow(ab_sel) == 0) {
    stop("No selected antibody-isoform rows found for ",
         target_cd, " / ", target_gene, " / ", target_ab)
  }
  
  tx_order <- ab_sel$transcript_id
  gene_seqs <- Biostrings::AAStringSet(ab_sel$protein_seq)
  names(gene_seqs) <- ab_sel$transcript_id
  
  if (length(gene_seqs) >= 2) {
    gene_msa <- msa::msa(gene_seqs, method = "ClustalOmega")
    gene_aligned <- msa_to_aastringset(gene_msa, expected_names = tx_order)
  } else {
    gene_aligned <- gene_seqs[tx_order]
  }
  
  tx_order <- names(gene_aligned)
  aln_len <- nchar(as.character(gene_aligned[[1]]))
  
  y_base <- setNames(rev(seq_along(tx_order)), tx_order)
  y_limits <- c(0.45, length(tx_order) + 0.75)
  
  tm_gene <- tm_seg2 %>%
    dplyr::mutate(transcript_id = strip_tx_version(transcript_id))
  
  if ("hgnc_symbol" %in% colnames(tm_gene)) {
    tm_gene <- tm_gene %>%
      dplyr::mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
      dplyr::filter(hgnc_symbol == target_gene)
  }
  
  tm_gene <- tm_gene %>%
    dplyr::filter(transcript_id %in% tx_order)
  
  aligned_topo <- lapply(tx_order, function(tid) {
    pep_len <- nchar(gsub("-", "", as.character(gene_aligned[[tid]])))
    seg_df <- tm_gene %>% dplyr::filter(transcript_id == tid)
    codes <- build_topology_codes(seg_df, pep_len)
    aln_seq <- as.character(gene_aligned[[tid]])
    aln_codes <- map_codes_to_alignment(aln_seq, codes)
    
    tibble::tibble(
      position = seq_along(aln_codes),
      transcript_id = tid,
      y = y_base[[tid]],
      topology_type = code_to_type(aln_codes)
    )
  }) %>%
    dplyr::bind_rows()
  
  backbone <- tibble::tibble(
    transcript_id = tx_order,
    y = unname(y_base[tx_order]),
    x0 = 1L,
    x1 = aln_len
  )
  
  required_pfam_cols <- c("transcript_id", "domain_label", "start", "end")
  missing_pfam_cols <- setdiff(required_pfam_cols, colnames(pfam_data))
  
  if (length(missing_pfam_cols) > 0) {
    warning(
      "Missing Pfam column(s): ",
      paste(missing_pfam_cols, collapse = ", "),
      ". Domains will not be drawn."
    )
    dom_raw <- tibble::tibble()
  } else {
    dom_raw <- pfam_data %>%
      dplyr::mutate(transcript_id = strip_tx_version(transcript_id)) %>%
      dplyr::filter(transcript_id %in% tx_order) %>%
      dplyr::transmute(
        transcript_id,
        domain_label = as.character(domain_label),
        aa_start = as.integer(start),
        aa_end = as.integer(end)
      ) %>%
      dplyr::filter(!is.na(aa_start), !is.na(aa_end), aa_start <= aa_end)
  }
  
  dom_aln <- project_domains_to_alignment(
    dom_raw = dom_raw,
    aligned = gene_aligned,
    y_base = y_base,
    ymin_off = -0.44,
    ymax_off = -0.25
  )
  
  prest_raw <- ab_sel %>%
    dplyr::filter(transcript_id %in% tx_order) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      interval = list(
        find_prest_interval(
          prest_seq = prest_seq,
          protein_seq = protein_seq,
          max_mismatch_frac = max_mismatch_frac
        )
      )
    ) %>%
    tidyr::unnest(interval) %>%
    dplyr::ungroup()
  
  interval_check <- prest_raw %>%
    dplyr::mutate(has_interval = !is.na(aa_start) & !is.na(aa_end)) %>%
    dplyr::count(present, has_interval, name = "n")
  
  message("\nInterval check for ", target_cd, " / ", target_gene, " / ", target_ab)
  print(interval_check)
  
  prest_aln <- project_prest_to_alignment(
    prest_raw = prest_raw,
    aligned = gene_aligned,
    y_base = y_base,
    ymin_off = 0.30,
    ymax_off = 0.50
  )
  
  topology_colors <- c(
    "Signal peptide" = "#C41C24",
    "Extracellular" = "#FFB20F",
    "Transmembrane" = "#18848C",
    "Intracellular" = "#96BDC6",
    "Alignment gap" = "#EDE7E3",
    "Unknown" = "grey70"
  )
  
  prest_color <- "#2D6A4F"
  structure_text_col <- "grey35"
  
  p_final <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = backbone,
      ggplot2::aes(x = x0, xend = x1, y = y, yend = y),
      linewidth = 0.35,
      color = "grey25"
    ) +
    ggplot2::geom_tile(
      data = aligned_topo,
      ggplot2::aes(x = position, y = y, fill = topology_type),
      height = 0.52,
      width = 1
    ) +
    ggplot2::scale_fill_manual(
      values = topology_colors,
      name = "Backbone topology (MSA bar)",
      drop = FALSE,
      guide = ggplot2::guide_legend(order = 2)
    ) +
    ggnewscale::new_scale_fill() +
    {
      if (nrow(dom_aln) > 0) {
        ggplot2::geom_rect(
          data = dom_aln,
          ggplot2::aes(
            xmin = aln_start,
            xmax = aln_end,
            ymin = ymin,
            ymax = ymax,
            fill = domain_label
          ),
          alpha = 0.9
        )
      } else {
        ggplot2::geom_blank()
      }
    } +
    {
      if (nrow(dom_aln) > 0) {
        ggplot2::scale_fill_brewer(
          palette = "Set2",
          name = "Pfam domain (below backbone)",
          labels = function(x) stringr::str_wrap(x, width = 45),
          guide = ggplot2::guide_legend(order = 3, ncol = 1)
        )
      } else {
        NULL
      }
    } +
    ggnewscale::new_scale_fill() +
    {
      if (nrow(prest_aln) > 0) {
        ggplot2::geom_rect(
          data = prest_aln,
          ggplot2::aes(
            xmin = aln_start,
            xmax = aln_end,
            ymin = ymin,
            ymax = ymax,
            fill = prest_status
          ),
          color = "grey25",
          linewidth = 0.15,
          alpha = 0.95
        )
      } else {
        ggplot2::geom_blank()
      }
    } +
    ggplot2::scale_fill_manual(
      values = c("HPA PrEST" = prest_color),
      name = "PrEST epitope (above backbone)",
      guide = ggplot2::guide_legend(order = 1)
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = seq(0, aln_len, by = 100)
    ) +
    ggplot2::scale_y_continuous(
      breaks = unname(y_base[tx_order]),
      labels = rep("", length(tx_order)),
      limits = y_limits,
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = "Amino-acid position in aligned protein isoforms",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = TEXT_SIZE) +
    ggplot2::theme(
      text = ggplot2::element_text(
        color = structure_text_col,
        family = "Helvetica",
        face = "plain",
        size = TEXT_SIZE
      ),
      axis.text.x = ggplot2::element_text(color = structure_text_col, size = TEXT_SIZE),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(color = structure_text_col, size = TEXT_SIZE),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "grey35", linewidth = 0.25),
      axis.ticks.x = ggplot2::element_line(color = "grey45", linewidth = 0.25),
      legend.title = ggplot2::element_text(color = structure_text_col, size = TEXT_SIZE, face = "plain"),
      legend.text = ggplot2::element_text(color = structure_text_col, size = TEXT_SIZE, face = "plain"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
  
  height_use <- max(2.7, 0.42 * length(tx_order) + 2)
  
  ggplot2::ggsave(
    output_file_pdf,
    p_final,
    width = 9.5,
    height = height_use,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  
  ggplot2::ggsave(
    sub("\\.pdf$", ".png", output_file_pdf),
    p_final,
    width = 9.5,
    height = height_use,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  
  invisible(list(
    plot = p_final,
    ab_sel = ab_sel,
    prest_raw = prest_raw,
    prest_aln = prest_aln,
    domain_aln = dom_aln,
    interval_check = interval_check
  ))
}

# ============================================================
# PANEL B: LOCALIZATION STRIP
# ============================================================

make_broad_loc_annotation <- function(loc_data,
                                      selected_tx,
                                      target_gene = NULL) {
  
  loc_tbl <- loc_data %>%
    dplyr::mutate(
      transcript_id = strip_tx_version(transcript_id),
      compartment = as.character(compartment),
      prob = as.numeric(prob),
      loc_class = dplyr::case_when(
        compartment == "Cell membrane" ~ "Membrane",
        compartment == "Extracellular" ~ "Secreted",
        TRUE ~ "Intracellular"
      )
    )
  
  if (!is.null(target_gene) && "hgnc_symbol" %in% colnames(loc_tbl)) {
    loc_tbl <- loc_tbl %>%
      dplyr::filter(hgnc_symbol == target_gene)
  }
  
  loc_tbl %>%
    dplyr::filter(
      transcript_id %in% selected_tx,
      !is.na(prob),
      !is.na(loc_class)
    ) %>%
    dplyr::group_by(transcript_id, loc_class) %>%
    dplyr::summarise(
      loc_prob = max(prob, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(
      dplyr::desc(loc_prob),
      factor(loc_class, levels = c("Membrane", "Secreted", "Intracellular"))
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::right_join(
      tibble::tibble(transcript_id = selected_tx),
      by = "transcript_id"
    ) %>%
    dplyr::mutate(
      loc_class = dplyr::coalesce(loc_class, "No prediction"),
      loc_class = factor(
        loc_class,
        levels = c("Membrane", "Secreted", "Intracellular", "No prediction")
      ),
      transcript_id = factor(transcript_id, levels = rev(selected_tx))
    )
}

make_locclass_strip <- function(loc_annot,
                                title = "Predicted\nlocalization") {
  
  plot_df <- loc_annot %>%
    dplyr::mutate(
      tx_label = as.character(transcript_id),
      loc_label = dplyr::case_when(
        loc_class == "No prediction" ~ "",
        TRUE ~ as.character(loc_class)
      )
    )
  
  ggplot2::ggplot(plot_df) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0.02, y = transcript_id, label = tx_label),
      hjust = 0,
      size = GEOM_TEXT_SIZE,
      color = text_col,
      family = "mono"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 2.05, y = transcript_id, label = loc_label),
      hjust = 0,
      size = GEOM_TEXT_SIZE,
      color = text_col,
      family = "sans"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 3.35),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(mult = c(0.14, 0.14))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    ggplot2::theme_void(base_size = TEXT_SIZE) +
    ggplot2::theme(
      text = ggplot2::element_text(size = TEXT_SIZE, color = text_col),
      plot.title = ggplot2::element_text(
        size = TITLE_SIZE,
        face = "bold",
        hjust = 0,
        color = text_col,
        lineheight = 0.95
      ),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

attach_loc_strip_to_structure <- function(loc_strip,
                                          structure_plot,
                                          panel_title = "PrEST and structure mapping") {
  
  structure_no_legend <- structure_plot +
    ggplot2::labs(title = panel_title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = TITLE_SIZE,
        face = "bold",
        hjust = 0.5,
        color = text_col
      ),
      legend.position = "none",
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
  
  structure_legend <- cowplot::get_legend(
    structure_plot +
      ggplot2::theme(
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        legend.key.height = grid::unit(0.42, "cm"),
        legend.key.width = grid::unit(0.42, "cm"),
        legend.title = ggplot2::element_text(size = TEXT_SIZE, color = text_col),
        legend.text = ggplot2::element_text(size = TEXT_SIZE, color = text_col)
      )
  )
  
  cowplot::plot_grid(
    NULL,
    loc_strip,
    structure_no_legend,
    structure_legend,
    ncol = 4,
    rel_widths = c(0.22, 1.75, 2.55, 1.95),
    align = "h",
    axis = "tb"
  )
}

# ============================================================
# PANEL C: FUNCTIONAL HEATMAP
# ============================================================

plot_tx_heatmap_top_union <- function(bundle_df,
                                      gene_symbol,
                                      population,
                                      selected_tx,
                                      tx_label_map,
                                      top_n = 8,
                                      padj_cut = 0.10,
                                      panel_title = "Transcript-associated biological programs",
                                      main_only = TRUE,
                                      drop_housekeeping = TRUE,
                                      extra_drop_pattern = NULL,
                                      select_positive_only = FALSE) {
  
  housekeeping_pattern <- paste(
    c(
      "small molecule",
      "amide metabolic process",
      "organic acid",
      "amino acid metabolic process",
      "carboxylic acid",
      "catabolic process",
      "metabolic process$",
      "actin filament organization",
      "ribonucleoprotein",
      "rna processing",
      "mrna processing",
      "protein folding",
      "proteasome",
      "ncRNA",
      "rRNA",
      "tRNA",
      "telomerase rna localization",
      "regulation of mitotic cytokinesis",
      "proton motive force",
      "atp synthesis coupled electron transport",
      "aerobic respiration",
      "oxidative phosphorylation",
      "cytoplasmic translation",
      "ribosome assembly",
      "ribosomal small subunit biogenesis"
    ),
    collapse = "|"
  )
  
  required_cols <- c(
    "gene", "population", "isoform", "pathway", "NES", "padj", "is_main_pathway"
  )
  
  missing_cols <- setdiff(required_cols, colnames(bundle_df))
  
  if (length(missing_cols) > 0) {
    stop(
      "function_main_bundle is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df <- bundle_df %>%
    dplyr::mutate(
      gene = as.character(gene),
      population = as.character(population),
      isoform = strip_tx_version(isoform),
      pathway_label = if ("pathway_pretty" %in% colnames(.)) {
        dplyr::coalesce(pathway_pretty, prettify_pathway(pathway))
      } else {
        prettify_pathway(pathway)
      }
    ) %>%
    dplyr::filter(
      .data$gene == .env$gene_symbol,
      .data$population == .env$population,
      .data$isoform %in% .env$selected_tx
    ) %>%
    dplyr::distinct(
      isoform,
      pathway,
      pathway_label,
      NES,
      padj,
      is_main_pathway,
      .keep_all = TRUE
    )
  
  if (nrow(df) == 0) {
    stop(
      "No rows found in function_main_bundle for gene = ",
      gene_symbol,
      ", population = ",
      population,
      "."
    )
  }
  
  if (main_only) {
    df <- df %>% dplyr::filter(is_main_pathway)
  }
  
  if (drop_housekeeping) {
    df <- df %>%
      dplyr::filter(
        !stringr::str_detect(
          stringr::str_to_lower(pathway_label),
          housekeeping_pattern
        )
      )
  }
  
  if (!is.null(extra_drop_pattern)) {
    df <- df %>%
      dplyr::filter(
        !stringr::str_detect(
          stringr::str_to_lower(pathway_label),
          stringr::str_to_lower(extra_drop_pattern)
        )
      )
  }
  
  df_for_selection <- df %>%
    dplyr::filter(!is.na(padj), padj <= padj_cut)
  
  if (select_positive_only) {
    df_for_selection <- df_for_selection %>%
      dplyr::filter(NES > 0)
  }
  
  top_union <- df_for_selection %>%
    dplyr::group_by(isoform) %>%
    dplyr::arrange(dplyr::desc(abs(NES)), padj, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(pathway, pathway_label)
  
  if (nrow(top_union) == 0) {
    stop(
      "No pathways found for ",
      gene_symbol,
      " after filtering. Try increasing padj_cut or setting main_only = FALSE."
    )
  }
  
  heat_df <- tidyr::expand_grid(
    isoform = selected_tx,
    pathway = unique(top_union$pathway)
  ) %>%
    dplyr::left_join(
      top_union %>% dplyr::distinct(pathway, pathway_label),
      by = "pathway"
    ) %>%
    dplyr::left_join(
      df %>% dplyr::select(isoform, pathway, pathway_label, NES, padj),
      by = c("isoform", "pathway", "pathway_label")
    ) %>%
    dplyr::mutate(
      isoform_label = dplyr::recode(
        isoform,
        !!!tx_label_map,
        .default = isoform
      )
    )
  
  pathway_order <- df_for_selection %>%
    dplyr::semi_join(top_union, by = c("pathway", "pathway_label")) %>%
    dplyr::group_by(pathway, pathway_label) %>%
    dplyr::summarise(
      max_abs_nes = max(abs(NES), na.rm = TRUE),
      dominant_tx = isoform[which.max(abs(NES))[1]],
      .groups = "drop"
    ) %>%
    dplyr::arrange(match(dominant_tx, selected_tx), dplyr::desc(max_abs_nes)) %>%
    dplyr::pull(pathway_label) %>%
    unique()
  
  heat_df <- heat_df %>%
    dplyr::mutate(
      pathway_label = factor(pathway_label, levels = rev(pathway_order)),
      isoform_label = factor(
        isoform_label,
        levels = unname(tx_label_map[selected_tx])
      )
    )
  
  nes_lim <- max(abs(heat_df$NES), na.rm = TRUE)
  if (!is.finite(nes_lim)) nes_lim <- 2
  
  ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = isoform_label, y = pathway_label)
  ) +
    ggplot2::geom_tile(
      fill = "grey92",
      color = "white",
      linewidth = 0.45
    ) +
    ggplot2::geom_tile(
      data = heat_df %>% dplyr::filter(!is.na(NES)),
      ggplot2::aes(fill = NES),
      color = "white",
      linewidth = 0.45
    ) +
    ggplot2::scale_fill_gradient2(
      low = nes_low,
      mid = nes_mid,
      high = nes_high,
      midpoint = 0,
      limits = c(-nes_lim, nes_lim),
      oob = scales::squish,
      na.value = "grey92",
      name = "NES"
    ) +
    ggplot2::labs(
      title = panel_title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = TEXT_SIZE) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      text = ggplot2::element_text(size = TEXT_SIZE, color = text_col),
      plot.title = ggplot2::element_text(
        size = TITLE_SIZE,
        face = "bold",
        hjust = 0.5,
        color = "grey20"
      ),
      axis.text.x = ggplot2::element_text(
        size = TEXT_SIZE,
        face = "plain",
        color = text_col,
        angle = 20,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = TEXT_SIZE, color = text_col),
      legend.title = ggplot2::element_text(size = TEXT_SIZE, color = text_col),
      legend.text = ggplot2::element_text(size = TEXT_SIZE, color = text_col),
      legend.position = "right",
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

# ============================================================
# SELECT FOCAL TRANSCRIPTS
# ============================================================

selected_tx_il7r_focus_top3 <- get_top_tx_for_gene_celltype(
  tx_summary = tx_summary,
  target_gene = "IL7R",
  focal_celltype = "CD4 T",
  top_n = 3
)

selected_tx_cd74_focus_top3 <- get_top_tx_for_gene_celltype(
  tx_summary = tx_summary,
  target_gene = "CD74",
  focal_celltype = "DC",
  top_n = 3
)

message("IL7R top 3 in CD4 T: ", paste(selected_tx_il7r_focus_top3, collapse = ", "))
message("CD74 top 3 in DC: ", paste(selected_tx_cd74_focus_top3, collapse = ", "))

tx_label_map_il7r_focus_top3 <- make_tx_label_map_for_selected(
  selected_tx = selected_tx_il7r_focus_top3,
  full_label_map = tx_label_map_il7r_all
)

tx_label_map_cd74_focus_top3 <- make_tx_label_map_for_selected(
  selected_tx = selected_tx_cd74_focus_top3,
  full_label_map = tx_label_map_cd74_all
)

# ============================================================
# FUNCTION BUNDLE SUBSETS
# ============================================================

function_main_bundle_IL7R_CD4T <- function_main_bundle %>%
  dplyr::filter(gene == "IL7R", population == "CD4_T")

function_main_bundle_CD74_DC <- function_main_bundle %>%
  dplyr::filter(gene == "CD74", population == "DC")

message("IL7R/CD4_T function rows: ", nrow(function_main_bundle_IL7R_CD4T))
message("CD74/DC function rows: ", nrow(function_main_bundle_CD74_DC))

# ============================================================
# PANEL A: TRANSCRIPT EXPRESSION LANDSCAPES
# ============================================================

pA_il7r <- make_transcript_landscape_top3_per_celltype(
  tx_summary = tx_summary,
  target_gene = "IL7R",
  gene_title = "",
  full_label_map = tx_label_map_il7r_all,
  top_n = 3,
  cell_order = cell_order,
  other_col = col_other
)

pA_cd74 <- make_transcript_landscape_top3_per_celltype(
  tx_summary = tx_summary,
  target_gene = "CD74",
  gene_title = "",
  full_label_map = tx_label_map_cd74_all,
  top_n = 3,
  cell_order = cell_order,
  other_col = col_other
)

# ============================================================
# PANEL B: STRUCTURE / LOCALIZATION / PrEST
# ============================================================

ab_tx_il7r_focus <- ab_tx %>%
  dplyr::filter(
    gene_symbol == "IL7R",
    ab_id == "HPA067550",
    transcript_id %in% selected_tx_il7r_focus_top3
  )

ab_tx_cd74_focus <- ab_tx %>%
  dplyr::filter(
    gene_symbol == "CD74",
    ab_id == "HPA010592",
    transcript_id %in% selected_tx_cd74_focus_top3
  )

res_il7r_structure <- plot_hpa_prest_structure(
  target_gene = "IL7R",
  target_cd = "CD127",
  target_ab = "HPA067550",
  ab_tx = ab_tx_il7r_focus,
  tm_seg2 = tm_seg2,
  pfam_data = ipr_df_pfam,
  output_file_pdf = file.path(
    structure_dir,
    "CD127_IL7R_HPA067550_PrEST_structure_top3_CD4T.pdf"
  ),
  selected_tx_keep = selected_tx_il7r_focus_top3,
  max_mismatch_frac = 0.05
)

res_cd74_structure <- plot_hpa_prest_structure(
  target_gene = "CD74",
  target_cd = "CD74",
  target_ab = "HPA010592",
  ab_tx = ab_tx_cd74_focus,
  tm_seg2 = tm_seg2,
  pfam_data = ipr_df_pfam,
  output_file_pdf = file.path(
    structure_dir,
    "CD74_HPA010592_PrEST_structure_top3_DC.pdf"
  ),
  selected_tx_keep = selected_tx_cd74_focus_top3,
  max_mismatch_frac = 0.05
)

loc_il7r_focus <- make_broad_loc_annotation(
  loc_data = cd_markers_localization,
  selected_tx = selected_tx_il7r_focus_top3,
  target_gene = "IL7R"
)

loc_cd74_focus <- make_broad_loc_annotation(
  loc_data = cd_markers_localization,
  selected_tx = selected_tx_cd74_focus_top3,
  target_gene = "CD74"
)

pB_il7r_loc_strip <- make_locclass_strip(
  loc_annot = loc_il7r_focus,
  title = "Predicted\nlocalization"
)

pB_cd74_loc_strip <- make_locclass_strip(
  loc_annot = loc_cd74_focus,
  title = "Predicted\nlocalization"
)

pB_il7r <- attach_loc_strip_to_structure(
  loc_strip = pB_il7r_loc_strip,
  structure_plot = res_il7r_structure$plot,
  panel_title = "CD127 / IL7R HPA067550 PrEST and protein-structure mapping"
)

pB_cd74 <- attach_loc_strip_to_structure(
  loc_strip = pB_cd74_loc_strip,
  structure_plot = res_cd74_structure$plot,
  panel_title = "PrEST and structure mapping"
)

# ============================================================
# PANEL C: FUNCTIONAL PROGRAM HEATMAPS
# ============================================================

pC_il7r <- plot_tx_heatmap_top_union(
  bundle_df = function_main_bundle_IL7R_CD4T,
  gene_symbol = "IL7R",
  population = "CD4_T",
  selected_tx = selected_tx_il7r_focus_top3,
  tx_label_map = tx_label_map_il7r_focus_top3,
  top_n = 8,
  padj_cut = 0.10,
  panel_title = "CD127 / IL7R transcript-associated programs in CD4 T",
  main_only = TRUE,
  drop_housekeeping = TRUE,
  extra_drop_pattern = NULL,
  select_positive_only = FALSE
)

pC_cd74 <- plot_tx_heatmap_top_union(
  bundle_df = function_main_bundle_CD74_DC,
  gene_symbol = "CD74",
  population = "DC",
  selected_tx = selected_tx_cd74_focus_top3,
  tx_label_map = tx_label_map_cd74_focus_top3,
  top_n = 8,
  padj_cut = 0.10,
  panel_title = "",
  main_only = TRUE,
  drop_housekeeping = TRUE,
  extra_drop_pattern = NULL,
  select_positive_only = FALSE
)

# ============================================================
# COMBINE PANELS
# ============================================================

fig_il7r <-
  patchwork::free(pA_il7r, side = "l") /
  patchwork::free(pB_il7r) /
  pC_il7r +
  patchwork::plot_layout(
    heights = c(0.90, 1.85, 1.35)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "CD127 / IL7R isoform landscape, protein structure, and functional programs"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = TAG_SIZE, color = "grey15"),
    plot.title = ggplot2::element_text(face = "bold", size = MAIN_TITLE_SIZE, hjust = 0, color = "grey15"),
    plot.margin = ggplot2::margin(3, 3, 3, 3)
  )

fig_cd74 <-
  patchwork::free(pA_cd74, side = "l") /
  patchwork::free(pB_cd74) /
  pC_cd74 +
  patchwork::plot_layout(
    heights = c(0.90, 2.05, 1.38)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = TAG_SIZE, color = "grey15"),
    plot.margin = ggplot2::margin(3, 3, 3, 3)
  )

fig_both <- fig_il7r | fig_cd74

# ============================================================
# SAVE GENE FIGURES
# ============================================================

save_plot_pair(
  fig_il7r,
  "IL7R_CD127_three_panel_workflow_top3_per_celltype_A_focal_top3_BC_locstrip",
  width = 12.2,
  height = 12.5,
  dpi = 600
)

save_plot_pair(
  fig_cd74,
  "CD74_three_panel_workflow_top3_per_celltype_A_focal_top3_BC_locstrip",
  width = 15.2,
  height = 13.5,
  dpi = 1200
)

save_plot_pair(
  fig_both,
  "IL7R_CD74_combined_three_panel_workflows_top3_per_celltype_A_focal_top3_BC_locstrip",
  width = 24,
  height = 12.8,
  dpi = 600
)

# ============================================================
# SAVE INDIVIDUAL PANELS
# ============================================================

save_plot_pair(
  pA_il7r,
  "PanelA_IL7R_transcript_landscape_top3_per_celltype",
  width = 8.8,
  height = 4.2,
  dpi = 600
)

save_plot_pair(
  pB_il7r,
  "PanelB_IL7R_PrEST_structure_focal_top3_CD4T_locstrip",
  width = 14.5,
  height = 4.8,
  dpi = 600
)

save_plot_pair(
  pC_il7r,
  "PanelC_IL7R_function_heatmap_CD4T_focal_top3",
  width = 8.8,
  height = 5.5,
  dpi = 600
)

save_plot_pair(
  pA_cd74,
  "PanelA_CD74_transcript_landscape_top3_per_celltype",
  width = 8.8,
  height = 4.2,
  dpi = 600
)

save_plot_pair(
  pB_cd74,
  "PanelB_CD74_PrEST_structure_focal_top3_DC_locstrip",
  width = 14.5,
  height = 4.8,
  dpi = 600
)

save_plot_pair(
  pC_cd74,
  "PanelC_CD74_function_heatmap_DC_focal_top3",
  width = 8.8,
  height = 5.5,
  dpi = 600
)

# ============================================================
# DISPLAY
# ============================================================

fig_il7r
fig_cd74
fig_both

message("Done. Outputs written to: ", fig_out_dir)
