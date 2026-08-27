# Tests for plot_utr_binding() and its helpers in R/plot_utr_binding.R

# --- .prep_utr_bed_tracks -------------------------------------------------

test_that(".prep_utr_bed_tracks reduces a single data frame and names it", {
  bed <- make_raw_bed(n = 2, chr = "chr1", start = c(100, 120),
                      end = c(150, 200), strand = "+")
  res <- .prep_utr_bed_tracks(bed, default_name = "mybed")
  expect_named(res, "mybed")
  expect_s4_class(res$mybed, "GRanges")
  expect_equal(length(res$mybed), 1L)                       # overlapping peaks merged
  expect_equal(GenomicRanges::start(res$mybed), 100)
  expect_equal(GenomicRanges::end(res$mybed), 200)
})

test_that(".prep_utr_bed_tracks keeps supplied names for a list of tracks", {
  res <- .prep_utr_bed_tracks(list(A = make_raw_bed(1), B = make_raw_bed(1, chr = "chr2")))
  expect_named(res, c("A", "B"))
})

test_that(".prep_utr_bed_tracks rejects duplicate names and bad elements", {
  expect_error(.prep_utr_bed_tracks(list(A = make_raw_bed(1), A = make_raw_bed(1))),
               class = "rnapeaks_error_invalid_arg")
  expect_error(.prep_utr_bed_tracks(list(make_raw_bed(1), 42)),
               class = "rnapeaks_error_invalid_arg")
  expect_error(.prep_utr_bed_tracks(42), class = "rnapeaks_error_invalid_arg")
  expect_error(.prep_utr_bed_tracks(list()), class = "rnapeaks_error_invalid_arg")
})

# --- plot_utr_side_map ----------------------------------------------------

one_track_side <- function() {
  d <- data.frame(position_in_region = 1:100,
                  frequency = c(rep(0, 20), rep(1, 10), rep(0, 70)),
                  group = "A", n_events = 5L, stringsAsFactors = FALSE)
  d$moving_avg <- d$frequency
  d
}

test_that("plot_utr_side_map returns a ggplot for a single track", {
  p <- plot_utr_side_map(one_track_side(), event_schema_utr, utr_style(),
                         "utr5", title = "5' UTR")
  expect_s3_class(p, "ggplot")
})

test_that("plot_utr_side_map aborts when a named palette misses a track", {
  d <- rbind(one_track_side(),
             transform(one_track_side(), group = "B"))
  style <- utr_style(palette = c(A = "red"))   # no entry for track B
  expect_error(plot_utr_side_map(d, event_schema_utr, style, "utr5"),
               class = "rnapeaks_error_invalid_arg")
})

# --- plot_utr_binding entry point -----------------------------------------

test_that("plot_utr_binding reports missing bed / bad species via the error boundary", {
  err <- expect_error(suppressMessages(plot_utr_binding()))
  expect_match(conditionMessage(err), "Failed to generate UTR binding plot")
  expect_s3_class(err$parent, "rnapeaks_error_invalid_arg")

  err2 <- expect_error(
    suppressMessages(plot_utr_binding(make_raw_bed(), species = "nope"))
  )
  expect_s3_class(err2$parent, "rnapeaks_error_invalid_arg")
})

test_that("plot_utr_binding returns 5' and 3' plot/data pairs", {
  skip_if_not_installed("rtracklayer")
  gtf <- write_utr_gtf(); on.exit(unlink(gtf), add = TRUE)
  bed <- make_raw_bed(n = 2, chr = "chr1", start = c(110, 410),
                      end = c(130, 430), strand = "+", protein = "SRSF1")
  res <- suppressMessages(plot_utr_binding(bed, gtf = gtf))
  expect_s3_class(res$utr5$plot, "ggplot")
  expect_s3_class(res$utr3$plot, "ggplot")
  expect_s3_class(res$utr5$data, "data.frame")
  expect_equal(nrow(res$utr5$data), event_schema_utr$n_bins)
})
