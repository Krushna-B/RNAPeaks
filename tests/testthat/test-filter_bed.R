# Tests for filter_bed() in R/bed-prepare.R
#
# Contract:
#   Restrict a validated BED to a chromosome / strand / coordinate window and
#   merge nearby peaks PER TARGET.
#   - Filtering is by OVERLAP: a peak is kept if it overlaps [start, end] at
#     all. Peaks that straddle a window edge are CLIPPED to the window bounds.
#   - chr is normalised, so "chr1" and "1" refer to the same chromosome.
#   - Empty result -> NULL.
#   - Merged output carries a `group_name` column (= source target).
#   (Track selection via `include` now happens upstream in check_bed().)

# --- argument validation --------------------------------------------------

test_that("start / end must be single, non-NA, parseable numbers with end >= start", {
  bed <- make_checked_bed()
  expect_error(filter_bed(bed, "1", NA, 1000, "+", 0), class = "rnapeaks_error_invalid_arg")
  expect_error(filter_bed(bed, "1", "abc", 1000, "+", 0), class = "rnapeaks_error_invalid_arg")
  expect_error(filter_bed(bed, "1", 0, NA, "+", 0), class = "rnapeaks_error_invalid_arg")
  expect_error(filter_bed(bed, "1", 500, 100, "+", 0), class = "rnapeaks_error_invalid_arg")
  expect_error(filter_bed(bed, "1", -1, 1000, "+", 0), class = "rnapeaks_error_invalid_arg")
})

test_that("numeric strings are accepted for start / end", {
  bed <- make_checked_bed()
  out <- filter_bed(bed, "1", "0", "1000", "+", 0)
  expect_equal(nrow(out), 3)
})

test_that("chr, strand and collapse are validated", {
  bed <- make_checked_bed()
  expect_error(filter_bed(bed, c("1", "2"), 0, 1000, "+", 0), class = "rnapeaks_error_invalid_arg")
  expect_error(filter_bed(bed, "1", 0, 1000, ".", 0), class = "rnapeaks_error_invalid_arg")
  expect_error(filter_bed(bed, "1", 0, 1000, "+", -1), class = "rnapeaks_error_invalid_arg")
  expect_error(filter_bed(bed, "1", 0, 1000, "+", NA), class = "rnapeaks_error_invalid_arg")
})

# --- window filtering (overlap + clip) ------------------------------------

test_that("peaks fully inside the window are kept and tagged with group_name", {
  out <- filter_bed(make_checked_bed(), "1", 0, 1000, "+", 0)
  expect_equal(nrow(out), 3)
  expect_true("group_name" %in% colnames(out))
  expect_equal(unique(out$group_name), "SRSF1")
})

test_that("peaks with no overlap of the window are dropped", {
  # peaks at 100-150, 300-350, 500-550; window 200-1000 excludes the first.
  out <- filter_bed(make_checked_bed(), "1", 200, 1000, "+", 0)
  expect_equal(nrow(out), 2)
})

test_that("peaks straddling a window edge are clipped to the window bounds", {
  # single peak 100-150; window 120-1000 clips its left edge to 120.
  out <- filter_bed(make_checked_bed(start = 100, end = 150), "1", 120, 1000, "+", 0)
  expect_equal(nrow(out), 1)
  expect_equal(out$start, 121)  # GRanges 1-based start of clipped 0-based 120
  expect_equal(out$end, 150)

  # single peak 500-550; window 0-520 clips its right edge to 520.
  out2 <- filter_bed(make_checked_bed(start = 500, end = 550), "1", 0, 520, "+", 0)
  expect_equal(nrow(out2), 1)
  expect_equal(out2$end, 520)
})

test_that("only the requested strand is retained", {
  bed <- rbind(
    make_checked_bed(start = c(100, 300), end = c(150, 350), strand = "+"),
    make_checked_bed(start = 700, end = 750, strand = "-")
  )
  out <- filter_bed(bed, "1", 0, 1000, "+", 0)
  expect_equal(nrow(out), 2)
})

test_that("chr is normalised so chr1 matches 1", {
  out <- filter_bed(make_checked_bed(), "chr1", 0, 1000, "+", 0)
  expect_equal(nrow(out), 3)
})

# --- empty result ---------------------------------------------------------

test_that("a window containing no peaks returns NULL", {
  out <- suppressMessages(filter_bed(make_checked_bed(), "1", 0, 10, "+", 0))
  expect_null(out)
})

# --- per-target merging ---------------------------------------------------

test_that("collapse merges nearby peaks of the SAME target", {
  bed <- make_checked_bed(start = c(100, 160), end = c(150, 200), target = "SRSF1")
  merged <- filter_bed(bed, "1", 0, 1000, "+", 1000)
  separate <- filter_bed(bed, "1", 0, 1000, "+", 0)
  expect_equal(nrow(merged), 1)
  expect_equal(nrow(separate), 2)
})

test_that("merging never crosses target boundaries", {
  bed <- rbind(
    make_checked_bed(start = 100, end = 150, target = "A"),
    make_checked_bed(start = 160, end = 200, target = "B")
  )
  out <- filter_bed(bed, "1", 0, 1000, "+", 1000)
  expect_equal(nrow(out), 2)
  expect_setequal(out$group_name, c("A", "B"))
})
