#' Internal: build per-gene UTR event table from a GTF
#'
#' Splits each protein-coding transcript's UTR rows into 5' and 3' sides
#' (relative to CDS), keeps only the exonic UTR pieces, and per gene picks
#' the transcript with the longest 5' UTR and the
#' transcript with the longest 3' UTR.
#'
#' @param gtf Normalized GTF data frame (see `get_GTF()`).
#' @param transcripts Optional character vector restricting which transcripts
#'   contribute. Accepts Ensembl transcript ids (`ENST…`), gene ids (`ENSG…`),
#'   or gene symbols; gene-level ids expand to all of that gene's transcripts.
#'   When `NULL`, all protein-coding transcripts contribute.
#'
#' @return List with:
#'   \describe{
#'     \item{events}{Data frame, one row per gene with `gene_id`,
#'       `gene_name`, `chr`, `strand`, `tx5`, `tx3`, `utr5_len`,
#'       `utr3_len`.}
#'     \item{utr5_pieces, utr3_pieces}{Data frames of exonic UTR pieces
#'       with `event_idx` (row index into `events`), `chr`, `start`,
#'       `end`, `strand`. Sorted by `(event_idx, start)`.}
#'   }
#'
#' @keywords internal
#' @noRd
build_utr_events <- function(gtf, transcripts = NULL) {
  if (!is.null(transcripts)) {
    if (!is.character(transcripts) || anyNA(transcripts) ||
        length(transcripts) == 0L) {
      abort_invalid_arg(
        "{.arg transcripts} must be a non-empty character vector with no NAs."
      )
    }
    tx_ids <- resolve_transcript_ids(gtf, transcripts)
    gtf <- gtf[gtf$transcript_id %in% tx_ids, , drop = FALSE]
  }

  # Restrict to protein-coding transcripts when biotype is present.
  if ("gene_biotype" %in% colnames(gtf)) {
    gtf <- gtf[is.na(gtf$gene_biotype) | gtf$gene_biotype == "protein_coding",
               , drop = FALSE]
  }

  keep_rows <- gtf$type %in% c("CDS", "UTR")
  feat <- gtf[keep_rows,
              c("seqnames", "start", "end", "strand", "type",
                "gene_id", "gene_name", "transcript_id"),
              drop = FALSE]
  if (nrow(feat) == 0L) {
    abort_not_found("GTF has no CDS or UTR rows after filtering.")
  }

  # Per-transcript CDS bounds.
  cds_idx  <- feat$type == "CDS"
  cds_tx   <- feat$transcript_id[cds_idx]
  cds_min  <- tapply(feat$start[cds_idx], cds_tx, min)
  cds_max  <- tapply(feat$end[cds_idx],   cds_tx, max)

  utr_idx <- feat$type == "UTR"
  utr_tx  <- feat$transcript_id[utr_idx]
  has_cds <- utr_tx %in% names(cds_min)
  if (!any(has_cds)) {
    abort_not_found("No UTR rows have a matching CDS in the same transcript.")
  }

  utr <- data.frame(
    transcript_id = utr_tx[has_cds],
    gene_id       = feat$gene_id[utr_idx][has_cds],
    gene_name     = feat$gene_name[utr_idx][has_cds],
    seqnames      = as.character(feat$seqnames[utr_idx])[has_cds],
    strand        = as.character(feat$strand[utr_idx])[has_cds],
    start         = feat$start[utr_idx][has_cds],
    end           = feat$end[utr_idx][has_cds],
    stringsAsFactors = FALSE
  )
  utr$cds_min <- cds_min[utr$transcript_id]
  utr$cds_max <- cds_max[utr$transcript_id]

  is_plus <- utr$strand == "+"
  upstream_of_cds   <- utr$end   <= utr$cds_min
  downstream_of_cds <- utr$start >= utr$cds_max
  utr$side <- NA_character_
  utr$side[(is_plus  & upstream_of_cds)   |
           (!is_plus & downstream_of_cds)] <- "utr5"
  utr$side[(is_plus  & downstream_of_cds) |
           (!is_plus & upstream_of_cds)]   <- "utr3"
  utr <- utr[!is.na(utr$side), , drop = FALSE]
  if (nrow(utr) == 0L) {
    abort_not_found("No UTR rows could be assigned to a 5' or 3' side.")
  }

  utr$width <- utr$end - utr$start + 1L

  # Per-(transcript, side) exonic length.
  is5 <- utr$side == "utr5"
  len5 <- tapply(utr$width[is5],  utr$transcript_id[is5],  sum)
  len3 <- tapply(utr$width[!is5], utr$transcript_id[!is5], sum)
  len_df <- data.frame(
    transcript_id = c(names(len5), names(len3)),
    side          = c(rep("utr5", length(len5)), rep("utr3", length(len3))),
    exonic_len    = as.integer(c(len5, len3)),
    stringsAsFactors = FALSE
  )

  # Transcript-level metadata.
  meta_idx  <- !duplicated(utr$transcript_id)
  tx_names  <- utr$transcript_id[meta_idx]
  meta_gene <- stats::setNames(utr$gene_id[meta_idx],   tx_names)
  meta_gn   <- stats::setNames(utr$gene_name[meta_idx], tx_names)
  meta_chr  <- stats::setNames(utr$seqnames[meta_idx],  tx_names)
  meta_str  <- stats::setNames(utr$strand[meta_idx],    tx_names)
  len_df$gene_id   <- meta_gene[len_df$transcript_id]
  len_df$gene_name <- meta_gn  [len_df$transcript_id]
  len_df$seqnames  <- meta_chr [len_df$transcript_id]
  len_df$strand    <- meta_str [len_df$transcript_id]

  # Pick the longest exonic transcript per (gene, side).
  len_df <- len_df[order(len_df$gene_id, len_df$side,
                         -len_df$exonic_len), , drop = FALSE]
  best <- len_df[!duplicated(len_df[, c("gene_id", "side")]), , drop = FALSE]

  best5 <- best[best$side == "utr5", , drop = FALSE]
  best3 <- best[best$side == "utr3", , drop = FALSE]

  gene_ids <- sort(unique(c(best5$gene_id, best3$gene_id)))
  i5 <- match(gene_ids, best5$gene_id)
  i3 <- match(gene_ids, best3$gene_id)

  base_meta <- rbind(best5[, c("gene_id", "gene_name", "seqnames", "strand")],
                     best3[, c("gene_id", "gene_name", "seqnames", "strand")])
  base_meta <- base_meta[!duplicated(base_meta$gene_id), , drop = FALSE]
  rownames(base_meta) <- base_meta$gene_id

  events <- data.frame(
    gene_id   = gene_ids,
    gene_name = base_meta[gene_ids, "gene_name"],
    chr       = base_meta[gene_ids, "seqnames"],
    strand    = base_meta[gene_ids, "strand"],
    tx5       = best5$transcript_id[i5],
    tx3       = best3$transcript_id[i3],
    utr5_len  = best5$exonic_len[i5],
    utr3_len  = best3$exonic_len[i3],
    stringsAsFactors = FALSE
  )
  events$utr5_len[is.na(events$utr5_len)] <- 0L
  events$utr3_len[is.na(events$utr3_len)] <- 0L

  # Flat per-side piece tables. event_idx is the row position in events.
  list(
    events      = events,
    utr5_pieces = .build_pieces(utr, "utr5", events, "tx5"),
    utr3_pieces = .build_pieces(utr, "utr3", events, "tx3")
  )
}

# Build a flat piece table for one side. Joins UTR rows to the events
# table via the chosen transcript per gene, attaches event_idx, and
# sorts by (event_idx, start) for per-event slicing downstream.
.build_pieces <- function(utr, side_label, events, tx_col) {
  picks <- utr[utr$side == side_label, , drop = FALSE]
  if (nrow(picks) == 0L) {
    return(data.frame(
      event_idx = integer(0),
      chr       = character(0),
      start     = integer(0),
      end       = integer(0),
      strand    = character(0),
      stringsAsFactors = FALSE
    ))
  }
  evt_idx <- match(picks$transcript_id, events[[tx_col]])
  keep    <- !is.na(evt_idx)
  out <- data.frame(
    event_idx = evt_idx[keep],
    chr       = picks$seqnames[keep],
    start     = picks$start[keep],
    end       = picks$end[keep],
    strand    = picks$strand[keep],
    stringsAsFactors = FALSE
  )
  out[order(out$event_idx, out$start), , drop = FALSE]
}


# Internal: resolve a mixed vector of user ids to GTF transcript_id values.
#
# Mirrors the id auto-detection in select_gene_transcript() so utr-binding
# accepts the same identifiers plot_gene() does. Each id is classified by shape:
# an Ensembl transcript id matches transcript_id directly; an Ensembl gene id or
# a gene symbol is expanded to *all* transcripts of that gene. Unmatched ids are
# reported but do not abort as long as at least one id resolves.
resolve_transcript_ids <- function(gtf, ids) {
  species  <- detect_species(gtf)
  tx_pat   <- if (species == "human") "^ENST\\d"  else "^ENSMUST\\d"
  gene_pat <- if (species == "human") "^ENSG\\d"  else "^ENSMUSG\\d"

  ids   <- trimws(ids)
  ids_u <- toupper(ids)
  kind  <- ifelse(grepl(tx_pat, ids_u), "tx",
                  ifelse(grepl(gene_pat, ids_u), "gene", "symbol"))

  out       <- character(0)
  unmatched <- character(0)
  for (i in seq_along(ids)) {
    tx <- switch(kind[i],
      tx     = gtf$transcript_id[!is.na(gtf$transcript_id) &
                                   gtf$transcript_id == ids_u[i]],
      gene   = gtf$transcript_id[!is.na(gtf$gene_id) &
                                   gtf$gene_id == ids_u[i]],
      symbol = gtf$transcript_id[!is.na(gtf$gene_name) &
                                   gtf$gene_name == normalize_label(ids[i], species)]
    )
    if (length(tx)) out <- c(out, tx) else unmatched <- c(unmatched, ids[i])
  }

  out <- unique(stats::na.omit(out))
  if (length(out) == 0L) {
    abort_not_found(c(
      "No GTF rows match {.arg transcripts}.",
      "i" = "Accepts Ensembl transcript ids ({.val ENST…}), gene ids ({.val ENSG…}), or gene symbols ({.val CXCR4})."
    ))
  }
  if (length(unmatched)) {
    cli::cli_alert_info("No match for: {.val {unmatched}}.")
  }
  out
}


# Gene-level sibling of resolve_transcript_ids(): tx / gene / symbol ids all
# collapse to gene_id. Returns unmatched instead of aborting so the caller can
# drop empty groups.
resolve_gene_ids <- function(gtf, ids) {
  species  <- detect_species(gtf)
  tx_pat   <- if (species == "human") "^ENST\\d"  else "^ENSMUST\\d"
  gene_pat <- if (species == "human") "^ENSG\\d"  else "^ENSMUSG\\d"

  ids   <- trimws(ids)
  ids_u <- toupper(ids)
  kind  <- ifelse(grepl(tx_pat, ids_u), "tx",
                  ifelse(grepl(gene_pat, ids_u), "gene", "symbol"))

  out       <- character(0)
  unmatched <- character(0)
  for (i in seq_along(ids)) {
    g <- switch(kind[i],
      tx     = gtf$gene_id[!is.na(gtf$transcript_id) &
                             gtf$transcript_id == ids_u[i]],
      gene   = gtf$gene_id[!is.na(gtf$gene_id) &
                             gtf$gene_id == ids_u[i]],
      symbol = gtf$gene_id[!is.na(gtf$gene_name) &
                             gtf$gene_name == normalize_label(ids[i], species)]
    )
    if (length(g)) out <- c(out, g) else unmatched <- c(unmatched, ids[i])
  }

  list(
    gene_ids  = unique(stats::na.omit(out)),
    unmatched = unmatched
  )
}
