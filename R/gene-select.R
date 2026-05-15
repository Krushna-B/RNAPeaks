#' Pull annotation rows for a single transcript from a GTF.
#'
#' Filters by `geneID` (if given), then picks the transcript: either an
#' explicit `TxID` or the longest transcript
#' for the gene. With only `TxID`, the gene filter is skipped and the
#' transcript is looked up directly in `gtf`.
#'
#' @param gtf GTF data frame.
#' @param geneID Gene symbol (any case) or Ensembl gene id.
#' @param TxID   Ensembl transcript id (any case). If `NULL`, the
#'   longest transcript for `geneID` is chosen.
#'
#' @return Data frame of GTF rows for the chosen transcript.
#' @noRd
#' @family gene
select_transcript <- function(gtf, geneID = NULL, TxID = NULL) {
  #Validate Params
  if (!is.data.frame(gtf) || nrow(gtf) == 0L) {
    abort_invalid_arg("{.arg gtf} must be a non-empty data frame.")
  }
  required_cols <- c("gene_id", "gene_name", "transcript_id", "type", "width")
  missing_cols  <- setdiff(required_cols, names(gtf))
  if (length(missing_cols)) {
    abort_invalid_gtf(
      "{.arg gtf} is missing required column{?s}: {.field {missing_cols}}"
    )
  }

  if (!is.null(geneID)) {
    if (length(geneID) != 1L) {
      abort_invalid_arg("{.arg geneID} must be length 1; got length {length(geneID)}.")
    }
    if (is.na(geneID) || !is.character(geneID) || !nzchar(geneID)) {
      abort_invalid_arg(
        "{.arg geneID} must be a single non-empty character string or {.code NULL}."
      )
    }
  }

  if (!is.null(TxID)) {
    if (length(TxID) != 1L) {
      abort_invalid_arg("{.arg TxID} must be length 1; got length {length(TxID)}.")
    }
    if (!is.na(TxID) && (!is.character(TxID) || !nzchar(TxID))) {
      abort_invalid_arg(
        "{.arg TxID} must be a single non-empty character string, {.code NA}, or {.code NULL}."
      )
    }
  }

  has_gene <- !is.null(geneID) && !is.na(geneID)
  has_tx   <- !is.null(TxID)   && !is.na(TxID)
  if (!has_gene && !has_tx) {
    abort_invalid_arg("Provide {.arg geneID} or {.arg TxID}.")
  }


  #Detect Species from the geneID col
  species <- detect_species(gtf)

  #Filter for gene id/name
  if (has_gene) {
    geneID_upper <- toupper(geneID)
    ens_pattern  <- if (species == "human") "^ENSG\\d" else "^ENSMUSG\\d"
    if (grepl(ens_pattern, geneID_upper)) {
      geneID_n <- geneID_upper
      field    <- "gene_id"
    } else {
      geneID_n <- normalize_label(geneID, species)
      field    <- "gene_name"
    }

    gtf <- gtf[!is.na(gtf[[field]]) & gtf[[field]] == geneID_n, , drop = FALSE]
    if (nrow(gtf) == 0L) {
      abort_not_found("Gene {.val {geneID}} not found in {species} GTF.")
    }
  }

  # Pick transcript
  if (has_tx) {
    tx_upper       <- toupper(TxID)
    ens_tx_pattern <- if (species == "human") "^ENST\\d" else "^ENSMUST\\d"

    if (grepl(ens_tx_pattern, tx_upper)) {
      # transcript id is uppercase for both species
      TxID_n   <- tx_upper
      tx_field <- "transcript_id"
    } else {
      # Transcript name is species specific casing
      if (!"transcript_name" %in% names(gtf)) {
        abort_invalid_gtf(c(
          "Cannot match by transcript name without {.field transcript_name} column.",
          "i" = "Pass an Ensembl transcript id ({.val ENST...} / {.val ENSMUST...}) instead."
        ))
      }
      TxID_n   <- normalize_label(TxID, species)
      tx_field <- "transcript_name"
    }

    if (!TxID_n %in% gtf[[tx_field]]) {
      in_gene <- if (has_gene) paste0(" in gene ", geneID) else ""
      abort_not_found("Transcript {.val {TxID}} not found{in_gene}.")
    }
  } else {
    txs <- unique(gtf[gtf$type == "transcript", c("transcript_id", "width")])
    if (nrow(txs) == 0L) {
      abort_not_found("No transcript rows for gene {.val {geneID}}.")
    }
    TxID_n   <- txs$transcript_id[which.max(txs$width)]
    tx_field <- "transcript_id"
  }

  # Filter to chosen transcript
  # Removes NA rows
  gtf[!is.na(gtf[[tx_field]]) & gtf[[tx_field]] == TxID_n, , drop = FALSE]
}


#' Pull annotation rows for every gene overlapping a genomic window.
#'
#' Keeps GTF rows for the longest transcript of each gene whose body
#' overlaps `chr:start-end` on `strand`.
#'
#' @param gtf GTF data frame.
#' @param chr Chromosome
#' @param start,end Window bounds in bp, `start <= end`.
#' @param strand `"+"` or `"-"`.
#'
#' @return Data frame of GTF rows for the selected transcripts.
#' @noRd
#' @family gene
select_region <- function(gtf, chr, start, end, strand) {
  # Validate gtf
  if (!is.data.frame(gtf) || nrow(gtf) == 0L) {
    abort_invalid_arg("{.arg gtf} must be a non-empty data frame.")
  }
  required_cols <- c("seqnames", "start", "end", "strand", "type",
                     "gene_id", "transcript_id", "width")
  missing_cols  <- setdiff(required_cols, names(gtf))
  if (length(missing_cols)) {
    abort_invalid_gtf(
      "{.arg gtf} is missing required column{?s}: {.field {missing_cols}}"
    )
  }

  # Validate window inputs
  check_string(chr, "chr")
  chr <- normalize_chr(chr)
  check_string(strand, "strand", choices = c("+", "-"))
  check_scalar_number(start, "start")
  check_scalar_number(end,   "end")
  if (start > end) {
    abort_invalid_arg(c(
      "{.arg start} must be <= {.arg end}.",
      "x" = "Got start = {.val {start}}, end = {.val {end}}."
    ))
  }

  # Filter to rows on the requested chromosome and strand overlapping the window
  seq_chr    <- as.character(gtf$seqnames)
  seq_strand <- as.character(gtf$strand)
  rows <- gtf[!is.na(seq_chr) & seq_chr == chr &
              !is.na(seq_strand) & seq_strand == strand &
              gtf$end >= start & gtf$start <= end, , drop = FALSE]
  if (nrow(rows) == 0L) {
    abort_not_found(
      "No features overlap {.val {chr}:{start}-{end}} on strand {.val {strand}}."
    )
  }

  # Restrict to transcript rows for the per-gene longest-transcript pick
  txs <- rows[rows$type == "transcript" & !is.na(rows$gene_id) &
              !is.na(rows$transcript_id) & !is.na(rows$width), , drop = FALSE]
  if (nrow(txs) == 0L) {
    abort_not_found(
      "No transcript rows found in {.val {chr}:{start}-{end}} on strand {.val {strand}}."
    )
  }

  # One transcript per gene: longest by width
  keep_ids <- vapply(
    split(seq_len(nrow(txs)), txs$gene_id),
    function(idx) txs$transcript_id[idx[which.max(txs$width[idx])]],
    character(1)
  )

  # Return all rows belonging to the selected transcripts
  rows[!is.na(rows$transcript_id) & rows$transcript_id %in% keep_ids,
       , drop = FALSE]
}


# Internal Helper:
#Detect which species from gtf
detect_species <- function(gtf) {
  ids <- stats::na.omit(gtf$gene_id)
  if (any(grepl("^ENSMUSG", ids))) return("mouse")
  if (any(grepl("^ENSG",    ids))) return("human")
  abort_species_unknown()
}

# Internal Helper:
# Normalize a label for species
#   human is all caps
#   mouse is title case
normalize_label <- function(label, species) {
  if (species == "human") {
    toupper(label)
  } else {
    paste0(toupper(substr(label, 1, 1)),
           tolower(substring(label, 2)))
  }
}
