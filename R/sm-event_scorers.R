#' Internal: score matrices for splicing / sequence map events
#'
#' Both scorers consume a regions `GRanges` from `schema$build_regions` and
#' return the same shape:
#' `[n_events, n_regions * region_width]`, each cell `0` or `1`.
#'
#' @keywords internal
#' @name event_scorers

# Place per-bp hits into the score matrix.
#
# `position_in_region` is always 1-based within its region. `position_in_transcript_order`
# tells us which orientation it's in:
#   FALSE (default): genomic order, e.g. from findOverlaps. Minus-strand positions
#     need to be mirrored within the region to land in plot/transcript order.
#   TRUE: already transcript order, e.g. motif starts from vmatchPattern on
#     getSeq output (getSeq already reverse-complemented minus-strand regions).
#     No within-region flip needed; only the region index is flipped.
.fill_score_matrix <- function(M, n_regions, region_width,
                               event_ids, region_indices,
                               strands, position_in_region,
                               position_in_transcript_order = FALSE) {
  is_minus <- strands == "-"

  plot_region_idx <- ifelse(is_minus,
                            n_regions - region_indices + 1L,
                            region_indices)
  plot_position_in_region <- if (position_in_transcript_order) {
    position_in_region
  } else {
    ifelse(is_minus,
           region_width - position_in_region + 1L,
           position_in_region)
  }

  col_idx <- (plot_region_idx - 1L) * region_width + plot_position_in_region

  in_bounds <- !is.na(col_idx)   & !is.na(event_ids) &
               col_idx   >= 1L   & col_idx   <= ncol(M) &
               event_ids >= 1L   & event_ids <= nrow(M)

  if (any(in_bounds)) {
    M[cbind(event_ids[in_bounds], col_idx[in_bounds])] <- 1L
  }
  M
}


#' @rdname event_scorers
#' @noRd
peaks_scorer <- function(regions_gr, bed_gr,
                         n_events, n_regions, region_width) {

  M <- matrix(0L, nrow = n_events, ncol = n_regions * region_width)
  if (length(regions_gr) == 0L || length(bed_gr) == 0L) return(M)

  overlaps <- GenomicRanges::findOverlaps(regions_gr, bed_gr,
                                          ignore.strand = FALSE)
  if (length(overlaps) == 0L) return(M)

  region_of_overlap <- S4Vectors::queryHits(overlaps)
  peak_of_overlap   <- S4Vectors::subjectHits(overlaps)

  intersections <- IRanges::pintersect(regions_gr[region_of_overlap],
                                       bed_gr[peak_of_overlap])
  widths <- IRanges::width(intersections)
  keep   <- widths > 0L
  if (!any(keep)) return(M)
  region_of_overlap <- region_of_overlap[keep]
  intersections     <- intersections[keep]
  widths            <- widths[keep]

  # Expand each intersection range to per-bp entries.
  intersection_of_bp <- rep(seq_along(intersections), widths)
  bp_genomic_pos     <- IRanges::start(intersections)[intersection_of_bp] +
                        sequence(widths) - 1L

  region_for_bp <- region_of_overlap[intersection_of_bp]
  event_ids     <- GenomicRanges::mcols(regions_gr)$event_id[region_for_bp]
  region_idxs   <- GenomicRanges::mcols(regions_gr)$region_idx[region_for_bp]
  strands       <- as.character(GenomicRanges::strand(regions_gr))[region_for_bp]
  region_starts <- GenomicRanges::start(regions_gr)[region_for_bp]

  .fill_score_matrix(
    M, n_regions, region_width,
    event_ids          = event_ids,
    region_indices     = region_idxs,
    strands            = strands,
    position_in_region = bp_genomic_pos - region_starts + 1L
  )
}


#' @rdname event_scorers
#' @noRd
motif_scorer <- function(regions_gr, genome, motifs,
                         n_events, n_regions, region_width) {

  M <- matrix(0L, nrow = n_events, ncol = n_regions * region_width)
  if (length(regions_gr) == 0L) return(M)

  regions_gr <- .align_seqnames_to_genome(regions_gr, genome)
  if (length(regions_gr) == 0L) return(M)

  region_seqs <- tryCatch(
    BSgenome::getSeq(genome, regions_gr),
    error = function(e) {
      cli::cli_warn(c(
        "Sequence extraction failed.",
        "x" = "{e$message}",
        "i" = "Check that {.arg genome} matches the chromosome naming of your events."
      ))
      NULL
    }
  )
  if (is.null(region_seqs)) return(M)

  region_widths  <- IRanges::width(regions_gr)
  region_strands <- as.character(GenomicRanges::strand(regions_gr))
  region_events  <- GenomicRanges::mcols(regions_gr)$event_id
  region_idxs    <- GenomicRanges::mcols(regions_gr)$region_idx

  for (motif in motifs) {
    hits <- Biostrings::vmatchPattern(
      Biostrings::DNAString(motif), region_seqs, fixed = FALSE
    )

    n_per_region <- lengths(hits)
    if (sum(n_per_region) == 0L) next

    region_of_hit <- rep(seq_along(regions_gr), n_per_region)
    starts        <- unlist(IRanges::start(hits), use.names = FALSE)

    # Drop hits past the actual region width.
    within <- starts <= region_widths[region_of_hit]
    if (!any(within)) next
    region_of_hit <- region_of_hit[within]
    starts        <- starts[within]

    M <- .fill_score_matrix(
      M, n_regions, region_width,
      event_ids                    = region_events[region_of_hit],
      region_indices               = region_idxs[region_of_hit],
      strands                      = region_strands[region_of_hit],
      position_in_region           = starts,
      position_in_transcript_order = TRUE
    )
  }
  M
}


# Map regions_gr chr names to the BSgenome's (e.g. "1" -> "chr1"), dropping
# anything the genome doesn't recognize.
.align_seqnames_to_genome <- function(regions_gr, genome) {
  genome_chrs <- GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(genome))
  levs        <- GenomeInfoDb::seqlevels(regions_gr)
  new_levs    <- ifelse(levs %in% genome_chrs, levs, paste0("chr", levs))

  unmapped <- !(new_levs %in% genome_chrs)
  if (any(unmapped)) {
    cli::cli_warn(
      "Dropping region{?s} on chromosome{?s} not in {.arg genome}: {.val {levs[unmapped]}}."
    )
  }

  regions_gr <- GenomeInfoDb::renameSeqlevels(
    regions_gr, stats::setNames(new_levs, levs)
  )
  GenomeInfoDb::seqlevels(regions_gr, pruning.mode = "coarse") <-
    intersect(GenomeInfoDb::seqlevels(regions_gr), genome_chrs)
  regions_gr
}
