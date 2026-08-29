# Tests for the exported sequence-map entry points:
#   skipped_exon_sequence_map, retained_intron_sequence_map,
#   five_prime_sequence_map (a5ss), three_prime_sequence_map (a3ss).
#
# These add motif normalization, genome resolution, and sequence extraction on
# top of the shared pipeline, so the happy paths need a BSgenome and are gated
# behind skip_if_not_installed(). The argument / motif validation runs before
# the genome is touched, so those error checks run unconditionally.
#
# Events use real hg38 chr1 coordinates so getSeq() succeeds; the motifs are not
# expected to occur there (the region is telomeric N's) -- the smoke tests check
# the returned structure, which holds regardless of hit count.

skip_no_hg38 <- function() skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

group_stats <- function() data.frame(
  GeneID = c("G1", "G2", "G3"),
  PValue = c(0.01, 0.01, 0.99),
  FDR    = c(0.01, 0.01, 0.99),
  IncLevelDifference = c(0.30, -0.30, 0.001),
  stringsAsFactors = FALSE
)

se_events <- function() cbind(
  data.frame(chr = "chr1", strand = "+"), group_stats(),
  data.frame(exonStart_0base = 14000, exonEnd = 15000,
             upstreamES = 10000, upstreamEE = 11000,
             downstreamES = 18000, downstreamEE = 19000)
)

ri_events <- function() cbind(
  data.frame(chr = "chr1", strand = "+"), group_stats(),
  data.frame(upstreamES = 10000, upstreamEE = 11000,
             downstreamES = 18000, downstreamEE = 19000)
)

# a5ss (+): alt exon left, flanking exon right.
a5_events <- function() cbind(
  data.frame(chr = "chr1", strand = "+"), group_stats(),
  data.frame(longExonStart_0base = 10000, longExonEnd = 12000,
             shortES = 10000, shortEE = 11000,
             flankingES = 15000, flankingEE = 16000)
)

# a3ss (+): flanking exon left, alt exon right.
a3_events <- function() cbind(
  data.frame(chr = "chr1", strand = "+"), group_stats(),
  data.frame(longExonStart_0base = 14000, longExonEnd = 16000,
             shortES = 15000, shortEE = 16000,
             flankingES = 10000, flankingEE = 11000)
)

seq_opts <- function() splicing_options(min_count = 0, width_exon = 10,
                                        width_intron = 20)

expect_map_result <- function(res) {
  expect_named(res, c("plot", "data"))
  expect_s3_class(res$plot, "ggplot")
  invisible(res)
}

# --- happy paths (need a genome) ------------------------------------------

test_that("skipped_exon_sequence_map (combined) returns one plot + data", {
  skip_no_hg38()
  res <- suppressMessages(
    skipped_exon_sequence_map(se_events(), "GCATG", opts = seq_opts())
  )
  expect_map_result(res)
})

test_that("individual mode returns one named result per motif", {
  skip_no_hg38()
  res <- suppressMessages(
    skipped_exon_sequence_map(se_events(), c("GCATG", "TGCAT"),
                              opts = seq_opts(), motif_mode = "individual")
  )
  expect_named(res, c("GCATG", "TGCAT"))
  expect_map_result(res$GCATG)
  expect_map_result(res$TGCAT)
})

test_that("retained_intron / 5' / 3' sequence maps run end to end", {
  skip_no_hg38()
  expect_map_result(suppressMessages(
    retained_intron_sequence_map(ri_events(), "GCATG", opts = seq_opts())))
  expect_map_result(suppressMessages(
    five_prime_sequence_map(a5_events(), "GCATG", opts = seq_opts())))
  expect_map_result(suppressMessages(
    three_prime_sequence_map(a3_events(), "GCATG", opts = seq_opts())))
})

# --- validation that runs before the genome is resolved -------------------

test_that("sequence maps report a missing events / sequence argument", {
  expect_error(skipped_exon_sequence_map(sequence = "GCATG"),
               class = "rnapeaks_error_invalid_arg")
  expect_error(skipped_exon_sequence_map(se_events()),
               class = "rnapeaks_error_invalid_arg")
})

test_that("an invalid IUPAC motif is rejected before scoring", {
  expect_error(skipped_exon_sequence_map(se_events(), "GCATX"),
               class = "rnapeaks_error_invalid_arg")
})

test_that("an unknown motif_mode is rejected", {
  expect_error(
    skipped_exon_sequence_map(se_events(), "GCATG", motif_mode = "sideways"),
    class = "rnapeaks_error_invalid_arg"
  )
})

test_that("the wrapper prefixes the failing map name onto the error", {
  expect_error(skipped_exon_sequence_map(sequence = "GCATG"),
               regexp = "skipped exon sequence map")
})
