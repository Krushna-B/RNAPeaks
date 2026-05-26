#' Internal: UTR plot scorer
#'
#' For each event with at least one overlapping peak:
#'   1. drop peaks whose width is not fully covered by the event's exonic
#'      UTR pieces (any base in an intronic gap disqualifies the peak);
#'   2. project retained peaks onto spliced (mRNA-order) coordinates;
#'   3. binarize per-bp and resample to `n_bins` via the within-bin mean.
#'
#'
#' @param pieces Data frame from [build_utr_events()] (one of
#'   `utr5_pieces` / `utr3_pieces`).
#' @param bed_gr Reduced `GRanges` of peaks (single track).
#' @param n_events Total number of events with any exonic UTR on this
#'   side.
#' @param n_bins Integer number of utr bins.
#'
#' @return `list(density = numeric(n_bins), n = integer)` — `n` is
#'   `n_events`.
#'
#' @keywords internal
#' @noRd
score_utr_metagene <- function(pieces, bed_gr, n_events, n_bins = 100L) {
  #Validate Params
  check_scalar_int(n_bins, "n_bins", min = 2L)
  if (!methods::is(bed_gr, "GRanges")) {
    abort_invalid_arg("{.arg bed_gr} must be a GRanges.")
  }

  if (n_events == 0L || nrow(pieces) == 0L) {
    return(list(density = rep(0, n_bins), n = as.integer(n_events)))
  }

  #Find Hits
  all_gr <- GenomicRanges::GRanges(
    seqnames = pieces$chr,
    ranges   = IRanges::IRanges(start = pieces$start, end = pieces$end),
    strand   = pieces$strand
  )
  hits <- GenomicRanges::findOverlaps(bed_gr, all_gr, ignore.strand = FALSE)
  if (length(hits) == 0L) {
    return(list(density = rep(0, n_bins), n = as.integer(n_events)))
  }

  peak_hit  <- S4Vectors::queryHits(hits)
  piece_hit <- S4Vectors::subjectHits(hits)
  evt_hit   <- pieces$event_idx[piece_hit]

  # Per-event slices, both for the pieces table and for the hit pairs.
  pieces_by_evt <- split(seq_len(nrow(pieces)), pieces$event_idx)
  hits_by_evt   <- split(seq_along(peak_hit),   evt_hit)

  # Cache peak coordinate vectors once outside the loop.
  peak_s <- as.integer(GenomicRanges::start(bed_gr))
  peak_e <- as.integer(GenomicRanges::end(bed_gr))
  peak_w <- as.integer(IRanges::width(bed_gr))

  density_sum <- numeric(n_bins)
  for (e_key in names(hits_by_evt)) {
    h <- hits_by_evt[[e_key]]
    pi <- pieces_by_evt[[e_key]]
    if (is.null(pi)) next

    vec <- .score_one_event(
      piece_s = pieces$start[pi],
      piece_e = pieces$end[pi],
      strand  = pieces$strand[pi][1L],
      peak_indices = peak_hit[h],
      peak_s = peak_s, peak_e = peak_e, peak_w = peak_w,
      n_bins = n_bins
    )
    density_sum <- density_sum + vec
  }

  list(
    density = density_sum / n_events,
    n       = as.integer(n_events)
  )
}

# Score a single event with vector inputs only
.score_one_event <- function(piece_s, piece_e, strand,
                             peak_indices,
                             peak_s, peak_e, peak_w, n_bins) {
  L <- sum(piece_e - piece_s + 1L)
  if (L == 0L) return(numeric(n_bins))

  # Drop peaks that span a UTR-internal intronic gap. Peaks that just
  # extend past the UTR edge (into CDS or upstream of TSS) are kept
  # only their UTR-overlapping bp end up in cov_bp below
  uniq_peaks <- unique(peak_indices)
  if (length(uniq_peaks) == 0L) return(numeric(n_bins))

  n_pieces <- length(piece_s)
  if (n_pieces > 1L) {
    gap_s <- piece_e[-n_pieces] + 1L
    gap_e <- piece_s[-1L]       - 1L
    pk_s_v <- peak_s[uniq_peaks]
    pk_e_v <- peak_e[uniq_peaks]
    has_gap <- logical(length(uniq_peaks))
    for (g in seq_along(gap_s)) {
      has_gap <- has_gap | (pk_s_v <= gap_e[g] & pk_e_v >= gap_s[g])
    }
    kept_peaks <- uniq_peaks[!has_gap]
  } else {
    kept_peaks <- uniq_peaks
  }
  if (length(kept_peaks) == 0L) return(numeric(n_bins))

  # Spliced coordinate offsets per piece, in mRNA order.
  piece_w <- piece_e - piece_s + 1L
  if (strand == "-") {
    ord <- rev(seq_along(piece_s))
  } else {
    ord <- seq_along(piece_s)
  }
  offsets <- c(0L, cumsum(piece_w[ord])[-length(ord)])

  # Binary spliced-position coverage.
  cov_bp <- integer(L)
  for (pk in kept_peaks) {
    ps <- peak_s[pk]; pe <- peak_e[pk]
    for (j in seq_along(ord)) {
      r  <- ord[j]
      rs <- piece_s[r]; re <- piece_e[r]
      lo <- max(ps, rs); hi <- min(pe, re)
      if (lo > hi) next
      local <- if (strand == "-") (re - hi + 1L):(re - lo + 1L)
               else                (lo - rs + 1L):(hi - rs + 1L)
      cov_bp[offsets[j] + local] <- 1L
    }
  }

  .bin_mean(cov_bp, n_bins)
}

# Resample length-L 0/1 vector to n_bins.
# - L >= n_bins: each bin's value is the within-bin mean (down-sample).
# - L <  n_bins: every bin gets the bp value at its nearest position
#   (up-sample).
.bin_mean <- function(x, n_bins) {
  L <- length(x)
  if (L == 0L) return(rep(0, n_bins))
  if (L >= n_bins) {
    bin <- as.integer(ceiling(seq_len(L) * n_bins / L))
    bin[bin > n_bins] <- n_bins
    sums   <- tabulate(bin[x > 0L], nbins = n_bins)
    counts <- tabulate(bin,         nbins = n_bins)
    out <- numeric(n_bins)
    nz <- counts > 0L
    out[nz] <- sums[nz] / counts[nz]
    return(out)
  }
  # Up-sample: bin k samples the bp at position ceil(k * L / n_bins).
  pos <- as.integer(ceiling(seq_len(n_bins) * L / n_bins))
  pos[pos < 1L] <- 1L
  pos[pos > L]  <- L
  as.numeric(x[pos])
}
