# Tests for plot_peaks_pipeline() in R/plot_peaks_pipeline.R
#
# The shared pipeline is exercised WITHOUT loading a real GTF: transcripts are
# supplied directly from the make_gtf() fixture. This checks the two behaviours
# that belong to the pipeline itself: (1) it renders a ggplot when peaks fall in
# the window, and (2) it returns NULL when no peaks survive filtering, including
# the gene-mode case where the window is DERIVED from the transcript bounds.

enst1_rows <- function() make_gtf()[make_gtf()$transcript_id == "ENST00000001", ]

test_that("gene mode renders a ggplot when peaks fall within the transcript window", {
  bed <- make_raw_bed(n = 1, chr = "chr1", start = 200, end = 250,
                      strand = "+", protein = "SRSF1")
  out <- suppressMessages(plot_peaks_pipeline(enst1_rows(), bed, is_region = FALSE))
  expect_s3_class(out, "ggplot")
})

test_that("gene mode derives its window from the transcript bounds", {
  # ENST1 spans 100-1100; a peak at 1200-1250 sits just past the derived window
  # and must be filtered out, yielding NULL rather than a plot.
  bed <- make_raw_bed(n = 1, chr = "chr1", start = 1200, end = 1250,
                      strand = "+", protein = "SRSF1")
  out <- suppressMessages(plot_peaks_pipeline(enst1_rows(), bed, is_region = FALSE))
  expect_null(out)
})

test_that("region mode renders a ggplot across multiple transcripts", {
  g  <- make_gtf()
  tx <- g[g$transcript_id %in% c("ENST00000001", "ENST00000003"), ]
  bed <- make_raw_bed(n = 2, chr = "chr1", start = c(200, 5100), end = c(250, 5150),
                      strand = "+", protein = "SRSF1")
  out <- suppressMessages(plot_peaks_pipeline(
    tx, bed, is_region = TRUE, window = list(start = 0, end = 6000)))
  expect_s3_class(out, "ggplot")
})

test_that("the pipeline returns NULL when no peaks remain after filtering", {
  bed <- make_raw_bed(n = 1, chr = "chr2", start = 200, end = 250,  # wrong chromosome
                      strand = "+", protein = "SRSF1")
  out <- suppressMessages(plot_peaks_pipeline(enst1_rows(), bed, is_region = FALSE))
  expect_null(out)
})
