# Tests for plot_gene() in R/plot_gene.R
#
# plot_gene wraps the pipeline in a tryCatch that re-aborts with the unified
# message "Failed to generate plot." The underlying RNAPeaks error is carried on
# the thrown condition's $parent, so error tests inspect the parent's class.

test_that("missing required arguments abort through the unified error boundary", {
  err <- expect_error(suppressMessages(plot_gene()))
  expect_match(conditionMessage(err), "Failed to generate plot")
  expect_s3_class(err$parent, "rnapeaks_error_invalid_arg")

  err2 <- expect_error(suppressMessages(plot_gene(bed = make_raw_bed())))  # gene missing
  expect_s3_class(err2$parent, "rnapeaks_error_invalid_arg")
})

test_that("an invalid species is reported before any annotation is loaded", {
  err <- expect_error(
    suppressMessages(plot_gene(make_raw_bed(), gene = "AAA", species = "nope"))
  )
  expect_s3_class(err$parent, "rnapeaks_error_invalid_arg")
})

test_that("an unknown gene surfaces as a not-found error", {
  skip_if_not_installed("rtracklayer")
  gtf <- write_min_gtf(); on.exit(unlink(gtf), add = TRUE)
  err <- expect_error(
    suppressMessages(plot_gene(make_raw_bed(), gene = "NOPE", gtf = gtf))
  )
  expect_s3_class(err$parent, "rnapeaks_error_not_found")
})

test_that("plot_gene returns a ggplot for a gene with overlapping peaks", {
  skip_if_not_installed("rtracklayer")
  gtf <- write_min_gtf(); on.exit(unlink(gtf), add = TRUE)
  bed <- make_raw_bed(n = 2, chr = "chr1", start = c(200, 900), end = c(250, 950),
                      strand = "+", protein = "SRSF1")
  g <- suppressMessages(plot_gene(bed, gene = "AAA", gtf = gtf))
  expect_s3_class(g, "ggplot")
})
