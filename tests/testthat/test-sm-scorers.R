# Tests for R/sm-event_scorers.R
#
# Started with the two position-mapping helpers, which decide where each hit
# lands in the plot's flattened [1 .. n_regions*region_width] coordinate:
#   col_idx = (plot_region - 1) * region_width + plot_position
# On the - strand the region ORDER is reversed (n_regions - r + 1). For peaks
# the within-region position is ALSO reversed (region_width - p + 1); for motifs
# (position_in_transcript_order = TRUE) it is kept as-is.

# --- .hits_to_col_idx: + strand -------------------------------------------

test_that("+ strand maps straight through: (region-1)*width + position", {
  expect_equal(.hits_to_col_idx(1, "+", 3, n_regions = 4, region_width = 10), 3)
  expect_equal(.hits_to_col_idx(2, "+", 1, n_regions = 4, region_width = 10), 11)
  expect_equal(.hits_to_col_idx(4, "+", 10, n_regions = 4, region_width = 10), 40)
})

# --- .hits_to_col_idx: - strand (peaks -> region AND position reversed) ----

test_that("- strand reverses both region order and within-region position for peaks", {
  # region 1 / pos 3 -> plot region 4, plot pos 8 -> (4-1)*10 + 8 = 38
  expect_equal(.hits_to_col_idx(1, "-", 3, n_regions = 4, region_width = 10), 38)
  # region 4 / pos 10 -> plot region 1, plot pos 1 -> 1 (mirror of the + case)
  expect_equal(.hits_to_col_idx(4, "-", 10, n_regions = 4, region_width = 10), 1)
})

# --- .hits_to_col_idx: - strand (motif -> region reversed, position kept) --

test_that("- strand in transcript order reverses the region but keeps the position", {
  # region 1 / pos 3 -> plot region 4, plot pos 3 -> (4-1)*10 + 3 = 33
  expect_equal(
    .hits_to_col_idx(1, "-", 3, n_regions = 4, region_width = 10,
                     position_in_transcript_order = TRUE),
    33
  )
})

test_that(".hits_to_col_idx aborts on unresolved '*' strand", {
  expect_error(.hits_to_col_idx(1, "*", 1, n_regions = 4, region_width = 10),
               regexp = "unknown|unresolved")
})

# --- invariant (independent of the exact formula) -------------------------
# Whatever the reorientation, mapping every (region, position) must be a
# bijection onto 1..n_positions: no collisions, nothing out of range, all
# positions reachable. A conceptually-wrong mapping breaks this even if its
# arithmetic is internally consistent.

test_that(".hits_to_col_idx is a permutation of 1..n_positions on both strands", {
  n_regions <- 4L; region_width <- 10L; n_pos <- n_regions * region_width
  rr <- rep(seq_len(n_regions), each = region_width)
  pp <- rep(seq_len(region_width), times = n_regions)

  for (strand in c("+", "-")) {
    cols_peaks <- .hits_to_col_idx(rr, rep(strand, n_pos), pp,
                                   n_regions, region_width)
    expect_setequal(cols_peaks, seq_len(n_pos))          # every position hit once
    expect_equal(length(unique(cols_peaks)), n_pos)      # no collisions

    cols_motif <- .hits_to_col_idx(rr, rep(strand, n_pos), pp,
                                   n_regions, region_width,
                                   position_in_transcript_order = TRUE)
    expect_setequal(cols_motif, seq_len(n_pos))
  }
})

# --- peaks_scorer ---------------------------------------------------------
# One event with 2 regions of width 10: region 1 at chr1:100-109 (region_idx 1),
# region 2 at chr1:200-209 (region_idx 2). Peak positions are known, so the
# expected col_idx is derived from the genomic offset, not the code.

two_region_gr <- function(strand = "+") {
  GenomicRanges::GRanges(
    "1", IRanges::IRanges(c(100, 200), c(109, 209)), strand = strand,
    event_id = c(1L, 1L), region_idx = c(1L, 2L)
  )
}

test_that("a peak covers exactly its overlapping bp at the right + strand positions", {
  peak <- GenomicRanges::GRanges("1", IRanges::IRanges(102, 104), strand = "+")
  h <- peaks_scorer(two_region_gr("+"), peak, n_regions = 2, region_width = 10)
  # genomic 102-104 in region 1 -> position_in_region 3,4,5 -> col 3,4,5
  expect_setequal(h$col_idx, c(3, 4, 5))
  expect_true(all(h$event_id == 1L))
})

test_that("a peak in region 2 is offset by one region width", {
  peak <- GenomicRanges::GRanges("1", IRanges::IRanges(205, 206), strand = "+")
  h <- peaks_scorer(two_region_gr("+"), peak, n_regions = 2, region_width = 10)
  # region 2, position 6,7 -> (2-1)*10 + c(6,7) = 16,17
  expect_setequal(h$col_idx, c(16, 17))
})

test_that("a peak covering an entire region lights exactly region_width positions", {
  peak <- GenomicRanges::GRanges("1", IRanges::IRanges(100, 109), strand = "+")
  h <- peaks_scorer(two_region_gr("+"), peak, n_regions = 2, region_width = 10)
  expect_setequal(h$col_idx, 1:10)      # all of region 1
  expect_equal(length(h$col_idx), 10)   # no duplicates
})

test_that("only the bp inside the region are counted (pintersect clips the peak)", {
  peak <- GenomicRanges::GRanges("1", IRanges::IRanges(105, 115), strand = "+")
  h <- peaks_scorer(two_region_gr("+"), peak, n_regions = 2, region_width = 10)
  expect_setequal(h$col_idx, 6:10)      # 105-109 only; 110-115 is past the region
})

test_that("- strand peaks are reoriented through the scorer", {
  peak <- GenomicRanges::GRanges("1", IRanges::IRanges(102, 104), strand = "-")
  h <- peaks_scorer(two_region_gr("-"), peak, n_regions = 2, region_width = 10)
  # region 1 / pos 3,4,5 on - strand -> plot region 2, positions 8,7,6 -> 16,17,18
  expect_setequal(h$col_idx, c(16, 17, 18))
})

test_that("peaks_scorer respects strand and empty inputs", {
  minus_peak <- GenomicRanges::GRanges("1", IRanges::IRanges(102, 104), strand = "-")
  h <- peaks_scorer(two_region_gr("+"), minus_peak, n_regions = 2, region_width = 10)
  expect_equal(length(h$col_idx), 0)    # opposite strands do not overlap

  empty <- GenomicRanges::GRanges()
  expect_equal(length(peaks_scorer(empty, minus_peak, 2, 10)$col_idx), 0)
})

# --- motif_scorer ---------------------------------------------------------
# Motif hit positions come from Biostrings::vmatchPattern; positions are kept in
# transcript order (position_in_transcript_order = TRUE), region order reversed
# on the - strand.

test_that("motif_scorer records each match start as a position in the region", {
  regions <- GenomicRanges::GRanges("1", IRanges::IRanges(100, 109), strand = "+",
                                    event_id = 1L, region_idx = 1L)
  seqs <- Biostrings::DNAStringSet("ACGTACGTAC")   # ACGT at positions 1 and 5
  h <- motif_scorer(regions, seqs, "ACGT", n_regions = 1, region_width = 10)
  expect_setequal(h$col_idx, c(1, 5))
  expect_true(all(h$event_id == 1L))
})

test_that("motif_scorer reverses region order on the - strand but keeps position", {
  regions <- GenomicRanges::GRanges(
    "1", IRanges::IRanges(c(100, 200), c(109, 209)), strand = "-",
    event_id = c(1L, 1L), region_idx = c(1L, 2L)
  )
  seqs <- Biostrings::DNAStringSet(c("ACGTAAAAAA",   # region 1: hit at start 1
                                     "AAAAAACGTA"))  # region 2: hit at start 6
  h <- motif_scorer(regions, seqs, "ACGT", n_regions = 2, region_width = 10)
  # region 1 (idx 1) -> plot region 2, pos 1 -> 11 ; region 2 (idx 2) -> plot region 1, pos 6 -> 6
  expect_setequal(h$col_idx, c(6, 11))
})

test_that("motif_scorer returns empty when nothing matches", {
  regions <- GenomicRanges::GRanges("1", IRanges::IRanges(100, 109), strand = "+",
                                    event_id = 1L, region_idx = 1L)
  seqs <- Biostrings::DNAStringSet("AAAAAAAAAA")
  h <- motif_scorer(regions, seqs, "ACGT", n_regions = 1, region_width = 10)
  expect_equal(length(h$col_idx), 0)
})

# --- .dedupe_hits ---------------------------------------------------------

test_that(".dedupe_hits keeps one row per (event, col) pair", {
  out <- .dedupe_hits(event_id = c(1, 1, 2), col_idx = c(5, 5, 5), n_positions = 40)
  expect_equal(out$event_id, c(1, 2))   # the duplicate (1, 5) collapses
  expect_equal(out$col_idx,  c(5, 5))
})

test_that(".dedupe_hits drops out-of-range, NA, and non-positive ids", {
  expect_equal(.dedupe_hits(c(1, 2), c(5, 50), n_positions = 40)$col_idx, 5)  # 50 > 40
  expect_equal(.dedupe_hits(c(1, 2), c(NA, 6), n_positions = 40)$col_idx, 6)
  expect_equal(.dedupe_hits(c(0, 1), c(5, 6), n_positions = 40)$event_id, 1)  # id 0 dropped
})

test_that(".dedupe_hits returns empty structure for empty input", {
  out <- .dedupe_hits(integer(0), integer(0), n_positions = 40)
  expect_equal(length(out$event_id), 0)
  expect_equal(length(out$col_idx), 0)
})
