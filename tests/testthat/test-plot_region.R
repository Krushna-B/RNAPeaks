# Tests for plot_region() in R/plot_region.R
#
# Same unified error boundary as plot_gene: the thrown error reads "Failed to
# generate plot." and carries the RNAPeaks error on $parent.

test_that("missing required arguments abort through the unified error boundary", {
  err <- expect_error(suppressMessages(plot_region()))
  expect_match(conditionMessage(err), "Failed to generate plot")
  expect_s3_class(err$parent, "rnapeaks_error_invalid_arg")
})

test_that("an invalid species is reported before any annotation is loaded", {
  err <- expect_error(
    suppressMessages(plot_region(make_raw_bed(), chr = "1", start = 0, end = 100,
                                 strand = "+", species = "nope"))
  )
  expect_s3_class(err$parent, "rnapeaks_error_invalid_arg")
})

test_that("a window overlapping no gene surfaces as a not-found error", {
  skip_if_not_installed("rtracklayer")
  gtf <- write_min_gtf(); on.exit(unlink(gtf), add = TRUE)
  err <- expect_error(
    suppressMessages(plot_region(make_raw_bed(), chr = "2", start = 0, end = 1000,
                                 strand = "+", gtf = gtf))
  )
  expect_s3_class(err$parent, "rnapeaks_error_not_found")
})

test_that("plot_region returns a ggplot for a window with overlapping peaks", {
  skip_if_not_installed("rtracklayer")
  gtf <- write_min_gtf(); on.exit(unlink(gtf), add = TRUE)
  bed <- make_raw_bed(n = 2, chr = "chr1", start = c(200, 900), end = c(250, 950),
                      strand = "+", protein = "SRSF1")
  g <- suppressMessages(plot_region(bed, chr = "1", start = 0, end = 2000,
                                    strand = "+", gtf = gtf))
  expect_s3_class(g, "ggplot")
})
