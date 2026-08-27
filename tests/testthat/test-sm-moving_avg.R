# Tests for apply_moving_average() in R/sm-moving_avg.R
#
# Windowed mean of width `window` centered with half = window %/% 2, i.e. each
# position i averages positions [i - half, i + (window-half-1)], CLAMPED to the
# region it belongs to so smoothing never bleeds across region boundaries.
# Every expected value below is computed by hand from that definition.

# --- pass-through cases ---------------------------------------------------

test_that("window NULL / <= 1 returns x unchanged", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(apply_moving_average(x, NULL, n_regions = 1, region_width = 5), x)
  expect_equal(apply_moving_average(x, 1,    n_regions = 1, region_width = 5), x)
  expect_equal(apply_moving_average(x, 0,    n_regions = 1, region_width = 5), x)
})

# --- odd window, single region --------------------------------------------

test_that("a width-3 window is a centered mean, clamped at the ends", {
  x <- c(1, 2, 3, 4, 5)
  # i=1: mean(1,2); i=2: mean(1,2,3); i=3: mean(2,3,4); i=4: mean(3,4,5); i=5: mean(4,5)
  expect_equal(apply_moving_average(x, 3, n_regions = 1, region_width = 5),
               c(1.5, 2, 3, 4, 4.5))
})

# --- even window ----------------------------------------------------------

test_that("a width-2 window averages [i-1, i]", {
  x <- c(1, 2, 3, 4)
  expect_equal(apply_moving_average(x, 2, n_regions = 1, region_width = 4),
               c(1, 1.5, 2.5, 3.5))
})

# --- region boundaries ----------------------------------------------------

test_that("smoothing is clamped within each region and never crosses the boundary", {
  x <- c(1, 2, 3, 10, 20, 30)   # region 1: 1..3, region 2: 10..30
  out <- apply_moving_average(x, 3, n_regions = 2, region_width = 3)
  # region 1 uses only 1,2,3; region 2 only 10,20,30
  expect_equal(out, c(1.5, 2, 2.5, 15, 20, 25))
})

# --- validation -----------------------------------------------------------

test_that("apply_moving_average validates x type, length, and window", {
  expect_error(apply_moving_average("a", 3, 1, 5), class = "rnapeaks_error_invalid_arg")
  expect_error(apply_moving_average(1:4, 3, n_regions = 1, region_width = 5),
               class = "rnapeaks_error_invalid_arg")   # length != n_regions*region_width
  expect_error(apply_moving_average(1:5, c(2, 3), 1, 5), class = "rnapeaks_error_invalid_arg")
  expect_error(apply_moving_average(1:5, NA, 1, 5), class = "rnapeaks_error_invalid_arg")
})
