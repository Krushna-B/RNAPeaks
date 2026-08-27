# Tests for R/gene-bam_tracks.R
#
# resolve_bam_names + the y-scaling arithmetic of prepare_bam_tracks are tested
# deterministically (the latter by mocking compute_bam_coverage, so no BAM or
# Rsamtools is required). compute_bam_coverage's validation is tested directly;
# its actual coverage extraction gets one Rsamtools-guarded integration test.

# --- resolve_bam_names ----------------------------------------------------

test_that("resolve_bam_names fills blank names from the file basename", {
  expect_equal(names(resolve_bam_names(c("dir/sample1.bam"))), "sample1")
  expect_equal(names(resolve_bam_names(c(Foo = "x.bam"))), "Foo")
  mixed <- resolve_bam_names(c(A = "a.bam", "path/b.bam"))
  expect_equal(names(mixed), c("A", "b"))
})

# --- prepare_bam_tracks: scaling math (mocked coverage) -------------------

test_that("prepare_bam_tracks returns NULL for no BAM files", {
  expect_null(prepare_bam_tracks(NULL, "1", 1, 5, base_y = 0, style = peaks_plot_style()))
})

test_that("prepare_bam_tracks scales coverage into the band using the data max", {
  testthat::local_mocked_bindings(
    compute_bam_coverage = function(bam_path, chr, start, end) {
      data.frame(pos = 1:5, coverage = c(0, 5, 10, 5, 0))
    }
  )
  style <- peaks_plot_style(bam_track_height = 2, bam_ylim = NULL)
  res <- prepare_bam_tracks(c(A = "a.bam"), "1", 1, 5, base_y = 10, style = style)

  # y_hi = 10 (data max), y_lo = 0; y_top = floor + cov/y_hi * band_h
  expect_equal(res$ribbons$y_top, c(10, 11, 12, 11, 10))
  expect_equal(unique(res$ribbons$y_bot), 10)
  expect_equal(res$labels$y, 11)             # midpoint of [10, 12]
  expect_equal(res$labels$track, "A")
  expect_equal(res$scales$y, c(10, 12))
  expect_equal(res$scales$text, c("0", "10"))
  expect_equal(res$total_height, 2)
})

test_that("prepare_bam_tracks honours a fixed bam_ylim", {
  testthat::local_mocked_bindings(
    compute_bam_coverage = function(bam_path, chr, start, end) {
      data.frame(pos = 1:5, coverage = c(0, 5, 10, 5, 0))
    }
  )
  style <- peaks_plot_style(bam_track_height = 2, bam_ylim = c(0, 20))
  res <- prepare_bam_tracks(c(A = "a.bam"), "1", 1, 5, base_y = 10, style = style)
  # scaled against 20, not the data max of 10
  expect_equal(res$ribbons$y_top, c(10, 10.5, 11, 10.5, 10))
  expect_equal(res$scales$text, c("0", "20"))
})

test_that("prepare_bam_tracks stacks multiple tracks on shared y-scale", {
  testthat::local_mocked_bindings(
    compute_bam_coverage = function(bam_path, chr, start, end) {
      if (bam_path == "a.bam") data.frame(pos = 1:3, coverage = c(0, 10, 0))
      else                     data.frame(pos = 1:3, coverage = c(0, 4, 0))
    }
  )
  style <- peaks_plot_style(bam_track_height = 2, bam_ylim = NULL)
  res <- prepare_bam_tracks(c(A = "a.bam", B = "b.bam"), "1", 1, 3, base_y = 0, style = style)

  expect_setequal(unique(res$ribbons$y_bot), c(0, 2))  # track B sits one band up
  expect_setequal(res$labels$y, c(1, 3))               # band midpoints
  expect_equal(res$total_height, 4)
})

# --- compute_bam_coverage: validation (no Rsamtools needed) ----------------

test_that("compute_bam_coverage validates the path and index before reading", {
  expect_error(compute_bam_coverage(123, "1", 1, 10), class = "rnapeaks_error_invalid_arg")
  expect_error(compute_bam_coverage(tempfile(fileext = ".bam"), "1", 1, 10),
               class = "rnapeaks_error_not_found")

  # file exists but no accompanying .bai index
  bam <- tempfile(fileext = ".bam")
  file.create(bam)
  on.exit(unlink(bam), add = TRUE)
  expect_error(compute_bam_coverage(bam, "1", 1, 10), class = "rnapeaks_error_not_found")
})

# --- compute_bam_coverage: real extraction (guarded) ----------------------

test_that("compute_bam_coverage returns per-base coverage over the region", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicAlignments")

  sam <- tempfile(fileext = ".sam")
  writeLines(c(
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:1\tLN:1000",
    "r1\t0\t1\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
    "r2\t0\t1\t105\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII"
  ), sam)
  bam <- Rsamtools::asBam(sam, destination = tempfile(), overwrite = TRUE)
  on.exit(unlink(c(sam, bam, paste0(bam, ".bai"))), add = TRUE)

  cov <- compute_bam_coverage(bam, "1", 100, 114)
  expect_equal(cov$pos, 100:114)
  # r1 covers 100-109, r2 covers 105-114 -> overlap 105-109 has depth 2
  expect_equal(cov$coverage,
               c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1))
})
