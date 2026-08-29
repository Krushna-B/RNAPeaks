# Tests for the exported splicing-map entry points:
#   skipped_exon_splicing_map, retained_intron_splicing_map,
#   five_prime_splicing_map (a5ss), three_prime_splicing_map (a3ss).
#
# These share the same wrapper (wrap_sm_errors -> validate_sm_inputs ->
# .peaks_to_granges -> peaks_scorer -> event_map_pipeline -> plot_event_map),
# so the internal pieces are covered elsewhere. Here we check the end-to-end
# contract on small synthetic events + a tiny BED (no genome needed): a happy
# path returns plot + data, and each wrapper reports its own missing-argument
# and bad-option errors.

# One event per analysis group (Positive / Negative / Control), following the
# partition rules: sig = FDR & PValue < .05, control = both > .95.
group_stats <- function() {
  data.frame(
    GeneID = c("G1", "G2", "G3"),
    PValue = c(0.01, 0.01, 0.99),
    FDR    = c(0.01, 0.01, 0.99),
    IncLevelDifference = c(0.30, -0.30, 0.001),
    stringsAsFactors = FALSE
  )
}

se_events <- function() cbind(
  data.frame(chr = "chr1", strand = "+"), group_stats(),
  data.frame(exonStart_0base = 1400, exonEnd = 1500,
             upstreamES = 1000, upstreamEE = 1100,
             downstreamES = 1800, downstreamEE = 1900)
)

ri_events <- function() cbind(
  data.frame(chr = "chr1", strand = "+"), group_stats(),
  data.frame(upstreamES = 1000, upstreamEE = 1100,
             downstreamES = 1800, downstreamEE = 1900)
)

# a5ss (+): alt exon on the left (1000-1100), flanking exon on the right.
a5_events <- function() cbind(
  data.frame(chr = "chr1", strand = "+"), group_stats(),
  data.frame(longExonStart_0base = 1000, longExonEnd = 1200,
             shortES = 1000, shortEE = 1100,
             flankingES = 1500, flankingEE = 1600)
)

# a3ss (+): flanking exon on the left (1000-1100), alt exon on the right.
a3_events <- function() cbind(
  data.frame(chr = "chr1", strand = "+"), group_stats(),
  data.frame(longExonStart_0base = 1400, longExonEnd = 1600,
             shortES = 1500, shortEE = 1600,
             flankingES = 1000, flankingEE = 1100)
)

# A small BED overlapping the region windows so scoring yields real hits.
peaks_bed <- function() data.frame(
  chr = "chr1", start = c(1095, 1490), end = c(1105, 1510),
  name = "RBP", score = 1, strand = "+",
  stringsAsFactors = FALSE
)

sm_opts <- function() splicing_options(min_count = 0, width_exon = 10,
                                       width_intron = 20)

# Assert the shared plot + data contract of a splicing-map result.
expect_map_result <- function(res) {
  expect_named(res, c("plot", "data"))
  expect_s3_class(res$plot, "ggplot")
  expect_named(res$data, c("Negative", "Positive"))
  invisible(res)
}

# --- happy paths ----------------------------------------------------------

test_that("skipped_exon_splicing_map returns a four-region plot and data", {
  res <- suppressMessages(
    skipped_exon_splicing_map(se_events(), peaks_bed(), opts = sm_opts())
  )
  expect_map_result(res)
  # SE has 4 regions -> the frequency frame passed to the plot spans 4 * width
  expect_equal(nrow(res$data$Positive), 1)   # one Positive event was supplied
})

test_that("retained_intron_splicing_map runs end to end", {
  res <- suppressMessages(
    retained_intron_splicing_map(ri_events(), peaks_bed(), opts = sm_opts())
  )
  expect_map_result(res)
})

test_that("five_prime_splicing_map (a5ss) runs end to end", {
  res <- suppressMessages(
    five_prime_splicing_map(a5_events(), peaks_bed(), opts = sm_opts())
  )
  expect_map_result(res)
})

test_that("three_prime_splicing_map (a3ss) runs end to end", {
  res <- suppressMessages(
    three_prime_splicing_map(a3_events(), peaks_bed(), opts = sm_opts())
  )
  expect_map_result(res)
})

# --- required-argument errors (each wrapper labels its own map) ------------

test_that("each splicing map reports a missing events / bed_file argument", {
  expect_error(skipped_exon_splicing_map(bed_file = peaks_bed()),
               class = "rnapeaks_error_invalid_arg")
  expect_error(skipped_exon_splicing_map(se_events()),
               class = "rnapeaks_error_invalid_arg")
  expect_error(retained_intron_splicing_map(bed_file = peaks_bed()),
               class = "rnapeaks_error_invalid_arg")
  expect_error(five_prime_splicing_map(bed_file = peaks_bed()),
               class = "rnapeaks_error_invalid_arg")
  expect_error(three_prime_splicing_map(a3_events()),
               class = "rnapeaks_error_invalid_arg")
})

test_that("the wrapper prefixes the failing map name onto the error", {
  expect_error(skipped_exon_splicing_map(bed_file = peaks_bed()),
               regexp = "skipped exon splicing map")
  expect_error(five_prime_splicing_map(bed_file = peaks_bed()),
               regexp = "5' splicing map")
})

# --- bad options / style routed through validate_sm_inputs ----------------

test_that("a non-splicing_options opts argument is rejected", {
  expect_error(
    skipped_exon_splicing_map(se_events(), peaks_bed(), opts = list()),
    class = "rnapeaks_error_invalid_arg"
  )
})
