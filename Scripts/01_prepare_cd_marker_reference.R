#!/usr/bin/env Rscript

# ============================================================
# 01_prepare_cd_marker_reference.R
#
# Build a CD marker transcript reference table from:
#   1. Ensembl release 115 GTF
#   2. curated CD marker gene list
#
# Output:
#   data/processed/CD_marker_transcript_reference_ensembl115.tsv
# ============================================================

source("scripts/00_setup.R")

message("Preparing CD marker transcript reference")
message("Using Ensembl release: ", ensembl_release)

# ============================================================
# CHECK INPUTS
# ============================================================

if (!file.exists(gtf_path)) {
  stop("Missing GTF file: ", gtf_path)
}

if (!file.exists(cd_marker_master_path)) {
  stop("Missing CD marker master list: ", cd_marker_master_path)
}

# ============================================================
# LOAD CD MARKER MASTER LIST
# ============================================================

cd_marker_master <- readr::read_csv(
  cd_marker_master_path,
  show_col_types = FALSE
) %>%
  dplyr::transmute(
    cd_name = stringr::str_trim(`CD Marker`),
    gene_symbol = stringr::str_trim(`Gene Symbol`)
  ) %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
  dplyr::distinct()

message("CD marker genes in master list: ", dplyr::n_distinct(cd_marker_master$gene_symbol))

# ============================================================
# IMPORT ENSEMBL GTF
# ============================================================

message("Importing GTF: ", gtf_path)

gtf <- rtracklayer::import(gtf_path)

# ============================================================
# TRANSCRIPT-LEVEL ANNOTATION
# ============================================================

tx_annot_all <- as.data.frame(gtf) %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::transmute(
    gene_id = normalize_ensg(gene_id),
    gene_symbol = gene_name,
    ensembl_transcript_id = normalize_enst(transcript_id),
    transcript_id_versioned = transcript_id,
    transcript_biotype = transcript_biotype,
    seqnames = as.character(seqnames),
    start = start,
    end = end,
    strand = as.character(strand)
  ) %>%
  dplyr::filter(
    !is.na(gene_id), gene_id != "",
    !is.na(gene_symbol), gene_symbol != "",
    !is.na(ensembl_transcript_id), ensembl_transcript_id != ""
  ) %>%
  dplyr::distinct()

message("Transcript annotation rows in GTF: ", nrow(tx_annot_all))
message("Annotated genes in GTF: ", dplyr::n_distinct(tx_annot_all$gene_symbol))

# ============================================================
# RESTRICT TO CD MARKER GENES
# ============================================================

cd_tx_reference <- tx_annot_all %>%
  dplyr::inner_join(cd_marker_master, by = "gene_symbol") %>%
  dplyr::mutate(
    ensembl_release = ensembl_release,
    coding_class = dplyr::case_when(
      transcript_biotype == "protein_coding" ~ "protein_coding",
      transcript_biotype == "protein_coding_LoF" ~ "protein_coding_LoF",
      transcript_biotype == "protein_coding_CDS_not_defined" ~ "CDS_not_defined",
      TRUE ~ "non_protein_coding"
    ),
    protein_coding_capable = transcript_biotype %in% c(
      "protein_coding",
      "protein_coding_LoF"
    )
  ) %>%
  dplyr::arrange(gene_symbol, ensembl_transcript_id)

message("CD marker genes found in GTF: ", dplyr::n_distinct(cd_tx_reference$gene_symbol))
message("CD marker transcripts found in GTF: ", dplyr::n_distinct(cd_tx_reference$ensembl_transcript_id))

missing_cd_genes <- cd_marker_master %>%
  dplyr::anti_join(
    cd_tx_reference %>% dplyr::distinct(gene_symbol),
    by = "gene_symbol"
  )

if (nrow(missing_cd_genes) > 0) {
  warning(
    "Some CD marker genes from the master list were not found in the GTF: ",
    paste(missing_cd_genes$gene_symbol, collapse = ", ")
  )
}

# ============================================================
# SAVE OUTPUTS
# ============================================================

out_reference <- file.path(
  proc_dir,
  "CD_marker_transcript_reference_ensembl115.tsv"
)

out_missing <- file.path(
  proc_dir,
  "CD_marker_genes_missing_from_ensembl115_gtf.tsv"
)

readr::write_tsv(cd_tx_reference, out_reference)
readr::write_tsv(missing_cd_genes, out_missing)

message("Saved CD marker transcript reference: ", out_reference)
message("Saved missing-gene table: ", out_missing)

# ============================================================
# QUICK SUMMARY TABLE
# ============================================================

summary_by_gene <- cd_tx_reference %>%
  dplyr::group_by(cd_name, gene_symbol) %>%
  dplyr::summarise(
    n_annotated_transcripts = dplyr::n_distinct(ensembl_transcript_id),
    n_protein_coding_capable = sum(protein_coding_capable),
    n_non_protein_coding = sum(coding_class == "non_protein_coding"),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_annotated_transcripts), gene_symbol)

out_summary <- file.path(
  proc_dir,
  "CD_marker_transcript_reference_summary_by_gene_ensembl115.tsv"
)

readr::write_tsv(summary_by_gene, out_summary)

message("Saved summary table: ", out_summary)
