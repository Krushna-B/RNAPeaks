#' Internal: hit-list scorers for splicing / sequence map events
#'
#' Both scorers consume a regions `GRanges` from `schema$build_regions` and
#' return `list(event_id, col_idx)` — parallel integer vectors of the
#' (event, position) pairs that scored 1 in the original dense layout.
#' `col_idx` is in `[1, n_regions * region_width]`.
#'
#' @keywords internal
#' @name event_scorers


.hits_to_col_idx <- function(region_indices, strands, position_in_region,
                             n_regions, region_width,
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
  (plot_region_idx - 1L) * region_width + plot_position_in_region
}


.dedupe_hits <- function(event_id, col_idx, n_positions) {
  if (length(event_id) == 0L) {
    return(list(event_id = integer(0), col_idx = integer(0)))
  }
  ok <- !is.na(col_idx)  & !is.na(event_id) &
        col_idx  >= 1L   & col_idx  <= n_positions &
        event_id >= 1L
  event_id <- event_id[ok]
  col_idx  <- col_idx[ok]
  if (length(event_id) == 0L) {
    return(list(event_id = integer(0), col_idx = integer(0)))
  }
  key  <- (as.numeric(event_id) - 1) * n_positions + col_idx
  keep <- !duplicated(key)
  list(event_id = event_id[keep], col_idx = col_idx[keep])
}


#' @rdname event_scorers
#' @noRd
peaks_scorer <- function(regions_gr, bed_gr, n_regions, region_width) {
  empty <- list(event_id = integer(0), col_idx = integer(0))
  if (length(regions_gr) == 0L || length(bed_gr) == 0L) return(empty)

  overlaps <- GenomicRanges::findOverlaps(regions_gr, bed_gr,
                                          ignore.strand = FALSE)
  if (length(overlaps) == 0L) return(empty)

  region_of_overlap <- S4Vectors::queryHits(overlaps)
  peak_of_overlap   <- S4Vectors::subjectHits(overlaps)

  intersections <- IRanges::pintersect(regions_gr[region_of_overlap],
                                       bed_gr[peak_of_overlap])
  widths <- IRanges::width(intersections)
  keep   <- widths > 0L
  if (!any(keep)) return(empty)
  region_of_overlap <- region_of_overlap[keep]
  intersections     <- intersections[keep]
  widths            <- widths[keep]

  intersection_of_bp <- rep(seq_along(intersections), widths)
  bp_genomic_pos     <- IRanges::start(intersections)[intersection_of_bp] +
                        sequence(widths) - 1L

  region_for_bp <- region_of_overlap[intersection_of_bp]
  event_id      <- GenomicRanges::mcols(regions_gr)$event_id[region_for_bp]
  region_idxs   <- GenomicRanges::mcols(regions_gr)$region_idx[region_for_bp]
  strands       <- as.character(GenomicRanges::strand(regions_gr))[region_for_bp]
  region_starts <- GenomicRanges::start(regions_gr)[region_for_bp]

  col_idx <- .hits_to_col_idx(
    region_indices     = region_idxs,
    strands            = strands,
    position_in_region = bp_genomic_pos - region_starts + 1L,
    n_regions          = n_regions,
    region_width       = region_width
  )

  .dedupe_hits(event_id, col_idx, n_regions * region_width)
}


#' @rdname event_scorers
#' @noRd
motif_scorer <- function(regions_gr, region_seqs, motifs,
                         n_regions, region_width) {
  empty <- list(event_id = integer(0), col_idx = integer(0))
  if (length(regions_gr) == 0L || length(region_seqs) == 0L ||
      length(motifs) == 0L) return(empty)

  region_widths  <- IRanges::width(regions_gr)
  region_strands <- as.character(GenomicRanges::strand(regions_gr))
  region_events  <- GenomicRanges::mcols(regions_gr)$event_id
  region_idxs    <- GenomicRanges::mcols(regions_gr)$region_idx

  region_indices_by_motif <- vector("list", length(motifs))
  hit_starts_by_motif     <- vector("list", length(motifs))

  for (i in seq_along(motifs)) {
    motif_hits <- Biostrings::vmatchPattern(
      Biostrings::DNAString(motifs[i]), region_seqs, fixed = FALSE
    )

    n_hits_per_region <- lengths(motif_hits)
    if (sum(n_hits_per_region) == 0L) next

    region_of_each_hit <- rep.int(seq_along(region_seqs), n_hits_per_region)
    hit_starts         <- unlist(IRanges::start(motif_hits), use.names = FALSE)

    keep <- hit_starts <= region_widths[region_of_each_hit]
    if (!any(keep)) next

    region_indices_by_motif[[i]] <- region_of_each_hit[keep]
    hit_starts_by_motif[[i]]     <- hit_starts[keep]
  }

  region_of_hit <- unlist(region_indices_by_motif, use.names = FALSE)
  hit_starts    <- unlist(hit_starts_by_motif,     use.names = FALSE)

  if (length(region_of_hit) == 0L) return(empty)

  col_idx <- .hits_to_col_idx(
    region_indices               = region_idxs[region_of_hit],
    strands                      = region_strands[region_of_hit],
    position_in_region           = hit_starts,
    n_regions                    = n_regions,
    region_width                 = region_width,
    position_in_transcript_order = TRUE
  )

  .dedupe_hits(region_events[region_of_hit], col_idx,
               n_regions * region_width)
}


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
