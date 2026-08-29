# Tests for R/utr-scorer.R -- the UTR density calculation core.
#
# Every expected vector below is derived by hand from the algorithm:
#   1. peaks that span an intronic gap between exonic pieces are dropped;
#   2. surviving peaks are projected onto SPLICED (mRNA-order) coordinates,
#      so on the - strand the high-coordinate end is position 0;
#   3. the 0/1 spliced vector is resampled to n_bins.
# Using n_bins == L (== 100) gives a 1:1 bin mapping, so bin index == bp index.

# --- .short_utr_resample --------------------------------------------------

test_that(".short_utr_resample keeps each bp in place and zero-pads the rest", {
  expect_equal(.short_utr_resample(c(1, 0, 1), 5), c(1, 0, 1, 0, 0))
})

# --- .bin_mean ------------------------------------------------------------

test_that(".bin_mean down-samples with the within-bin mean when L >= n_bins", {
  expect_equal(.bin_mean(c(1, 1, 0, 0), 2), c(1, 0))
  expect_equal(.bin_mean(c(1, 0, 1, 0), 2), c(0.5, 0.5))
  expect_equal(.bin_mean(c(0, 0, 1, 1), 2), c(0, 1))
})

test_that(".bin_mean falls back to bp-in-place when L < n_bins, and handles L == 0", {
  expect_equal(.bin_mean(c(1, 1), 5), c(1, 1, 0, 0, 0))
  expect_equal(.bin_mean(integer(0), 3), c(0, 0, 0))
})

# --- .score_one_event: single-piece projection ----------------------------

test_that("a peak inside a single-piece + UTR lands on the matching bins", {
  vec <- .score_one_event(piece_s = 100, piece_e = 199, strand = "+",
                          peak_indices = 1,
                          peak_s = 120, peak_e = 129, peak_w = 10, n_bins = 100)
  expect_equal(which(vec == 1), 21:30)   # (120-100+1)..(129-100+1)
  expect_equal(sum(vec), 10)
})

test_that("- strand projection mirrors position (high coordinate = spliced 0)", {
  vec <- .score_one_event(piece_s = 100, piece_e = 199, strand = "-",
                          peak_indices = 1,
                          peak_s = 120, peak_e = 129, peak_w = 10, n_bins = 100)
  expect_equal(which(vec == 1), 71:80)   # (199-129+1)..(199-120+1)
})

# --- .score_one_event: multi-piece splicing + gap dropping -----------------

test_that("peaks are placed in spliced coordinates across two pieces", {
  # pieces 100-149 (offset 0) and 200-249 (offset 50), gap 150-199.
  in_p1 <- .score_one_event(c(100, 200), c(149, 249), "+", 1,
                            peak_s = 110, peak_e = 119, peak_w = 10, n_bins = 100)
  expect_equal(which(in_p1 == 1), 11:20)

  in_p2 <- .score_one_event(c(100, 200), c(149, 249), "+", 1,
                            peak_s = 210, peak_e = 219, peak_w = 10, n_bins = 100)
  expect_equal(which(in_p2 == 1), 61:70)   # offset 50 + (210-200+1)..(219-200+1)
})

test_that("a peak spanning the intronic gap between pieces is dropped", {
  vec <- .score_one_event(c(100, 200), c(149, 249), "+", 1,
                          peak_s = 140, peak_e = 210, peak_w = 71, n_bins = 100)
  expect_equal(sum(vec), 0)
})

# --- score_utr_side: overlap + averaging -----------------------------------

single_event_pieces <- function() {
  data.frame(event_idx = 1, chr = "1", start = 100, end = 199, strand = "+",
             stringsAsFactors = FALSE)
}
peak_gr <- function(s, e, strand = "+") {
  GenomicRanges::GRanges("1", IRanges::IRanges(s, e), strand = strand)
}

test_that("score_utr_side projects an overlapping peak and reports n_events", {
  s <- score_utr_side(single_event_pieces(), peak_gr(120, 129),
                      n_events = 1, n_bins = 100)
  expect_equal(which(s$density == 1), 21:30)
  expect_equal(s$n, 1L)
})

test_that("score_utr_side averages the summed signal over ALL events", {
  pieces <- rbind(
    data.frame(event_idx = 1, chr = "1", start = 100, end = 199, strand = "+"),
    data.frame(event_idx = 2, chr = "1", start = 300, end = 399, strand = "+")
  )
  # only event 1 has a peak, but the mean divides by n_events = 2
  s <- score_utr_side(pieces, peak_gr(120, 129), n_events = 2, n_bins = 100)
  expect_equal(unique(s$density[21:30]), 0.5)
  expect_equal(s$n, 2L)
})

test_that("score_utr_side respects strand when overlapping", {
  s <- score_utr_side(single_event_pieces(), peak_gr(120, 129, strand = "-"),
                      n_events = 1, n_bins = 100)
  expect_equal(sum(s$density), 0)   # - peak does not overlap the + piece
})

test_that("score_utr_side returns a flat zero curve for empty inputs", {
  expect_equal(score_utr_side(single_event_pieces(), peak_gr(120, 129),
                              n_events = 0, n_bins = 100)$density,
               rep(0, 100))
  empty_pieces <- single_event_pieces()[0, ]
  expect_equal(sum(score_utr_side(empty_pieces, peak_gr(120, 129),
                                  n_events = 3, n_bins = 100)$density), 0)
})

test_that("score_utr_side validates n_bins and bed_gr", {
  expect_error(score_utr_side(single_event_pieces(), peak_gr(1, 2), 1, n_bins = 1),
               class = "rnapeaks_error_invalid_arg")
  expect_error(score_utr_side(single_event_pieces(), "notgr", 1, n_bins = 100),
               class = "rnapeaks_error_invalid_arg")
})

# --- .score_one_side ------------------------------------------------------

test_that(".score_one_side builds a per-bin frame per (track, group)", {
  events <- data.frame(gene_id = "g1", utr5_len = 100L, stringsAsFactors = FALSE)
  df <- .score_one_side(list(A = peak_gr(120, 129)), single_event_pieces(),
                        events, group_events = list(`All genes` = 1L),
                        len_col = "utr5_len", n_bins = 100)
  expect_equal(nrow(df), 100)
  expect_equal(unique(df$track), "A")
  expect_equal(unique(df$gene_group), "All genes")
  expect_equal(which(df$frequency == 1), 21:30)
  expect_true(all(df$n_events == 1))
})

test_that(".score_one_side yields a row block per group x track", {
  events <- data.frame(gene_id = "g1", utr5_len = 100L, stringsAsFactors = FALSE)
  df <- .score_one_side(list(A = peak_gr(120, 129)), single_event_pieces(),
                        events,
                        group_events = list(High = 1L, Low = 1L),
                        len_col = "utr5_len", n_bins = 100)
  expect_equal(nrow(df), 200)
  expect_setequal(unique(df$gene_group), c("High", "Low"))
})

# --- .smooth_side ---------------------------------------------------------

test_that(".smooth_side is a no-op when smoothing is disabled", {
  df <- data.frame(track = "A", gene_group = "All genes",
                   frequency = c(0, 1, 0, 1, 0))
  expect_equal(.smooth_side(df, window = 0, n_bins = 5), df$frequency)
  expect_equal(.smooth_side(df, window = NULL, n_bins = 5), df$frequency)
  expect_equal(.smooth_side(df, window = 1, n_bins = 5), df$frequency)
})

test_that(".smooth_side preserves a constant curve", {
  df <- data.frame(track = "A", gene_group = "All genes",
                   frequency = rep(0.5, 20))
  expect_true(all(abs(.smooth_side(df, window = 3, n_bins = 20) - 0.5) < 1e-9))
})
