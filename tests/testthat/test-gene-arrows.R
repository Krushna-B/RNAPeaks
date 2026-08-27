# Tests for R/gene-arrows.R
#
# Verifies the arrow-placement arithmetic exactly:
#   n arrows -> spacing = vis_len / (n+1); centers = vis_start + spacing*k;
#   shaft = shaft_frac * spacing; a "+" strand arrow runs center-shaft/2 ..
#   center+shaft/2, a "-" strand arrow runs the other way.

# --- placement math (+ strand) --------------------------------------------

test_that("arrows are evenly spaced with correct shaft length on + strand", {
  introns <- data.frame(start = 0, end = 900, strand = "+", y = 0)
  out <- build_intron_arrows(introns, xlim = c(0, 900),
                             arrows_per_view = 12, max_per_intron = 2, shaft_frac = 0.4)
  # vis_len 900, n = min(2, round(1*12)) = 2 -> spacing 300, centers 300/600,
  # shaft = 0.4*300 = 120 -> +/- 60 around each center.
  expect_equal(nrow(out), 2)
  expect_equal(out$x,    c(240, 540))
  expect_equal(out$xend, c(360, 660))
  expect_equal(out$y,    c(0, 0))
})

test_that("- strand arrows point the opposite way", {
  introns <- data.frame(start = 0, end = 900, strand = "-", y = 0)
  out <- build_intron_arrows(introns, xlim = c(0, 900),
                             arrows_per_view = 12, max_per_intron = 2, shaft_frac = 0.4)
  expect_equal(out$x,    c(360, 660))   # center + shaft/2
  expect_equal(out$xend, c(240, 540))   # center - shaft/2
})

# --- density controls -----------------------------------------------------

test_that("arrow count is capped by max_per_intron", {
  introns <- data.frame(start = 0, end = 1000, strand = "+", y = 0)
  out <- build_intron_arrows(introns, xlim = c(0, 1000), arrows_per_view = 100)
  expect_equal(nrow(out), 6)   # default max_per_intron = 6
})

test_that("a visible-but-small intron still gets at least one arrow", {
  introns <- data.frame(start = 0, end = 30, strand = "+", y = 0)
  out <- build_intron_arrows(introns, xlim = c(0, 1000))
  # vis_frac 0.03 >= 0.02 threshold; round(0.03*12)=0 -> floored up to 1
  expect_equal(nrow(out), 1)
})

test_that("introns below min_intron_frac or min_intron_bp receive no arrows", {
  tiny_frac <- data.frame(start = 0, end = 5, strand = "+", y = 0)   # frac 0.005 < 0.02
  expect_equal(nrow(build_intron_arrows(tiny_frac, xlim = c(0, 1000))), 0)

  short_bp <- data.frame(start = 0, end = 50, strand = "+", y = 0)
  expect_equal(nrow(build_intron_arrows(short_bp, xlim = c(0, 1000), min_intron_bp = 100)), 0)
})

test_that("arrows are confined to the visible window", {
  introns <- data.frame(start = 0, end = 2000, strand = "+", y = 0)
  out <- build_intron_arrows(introns, xlim = c(500, 1500), max_per_intron = 2)
  centers <- (out$x + out$xend) / 2
  expect_equal(nrow(out), 2)
  expect_true(all(centers >= 500 & centers <= 1500))
})

# --- degenerate inputs ----------------------------------------------------

test_that("empty introns or a zero-width view yield an empty arrow frame", {
  empty <- data.frame(start = numeric(0), end = numeric(0),
                      strand = character(0), y = numeric(0))
  out <- build_intron_arrows(empty, xlim = c(0, 1000))
  expect_equal(nrow(out), 0)
  expect_equal(names(out), c("x", "xend", "y", "yend"))

  introns <- data.frame(start = 0, end = 900, strand = "+", y = 0)
  expect_equal(nrow(build_intron_arrows(introns, xlim = c(500, 500))), 0)  # view_span 0
})

# --- empty_arrows / resolve_intron_y --------------------------------------

test_that("empty_arrows has the canonical numeric columns and no rows", {
  ea <- empty_arrows()
  expect_equal(names(ea), c("x", "xend", "y", "yend"))
  expect_equal(nrow(ea), 0)
  expect_type(ea$x, "double")
})

test_that("resolve_intron_y prefers y, then y_mid, then y_start/y_end mean, else 0", {
  expect_equal(resolve_intron_y(data.frame(start = 0, end = 1, y = 5)), 5)
  expect_equal(resolve_intron_y(data.frame(start = 0, end = 1, y_mid = 3)), 3)
  expect_equal(resolve_intron_y(data.frame(y_start = 0, y_end = 1)), 0.5)
  expect_equal(resolve_intron_y(data.frame(start = 0, end = 1)), 0)
})

# --- validators -----------------------------------------------------------

test_that("require_intron_df / require_xlim / require_unit_fraction reject bad input", {
  expect_error(build_intron_arrows(list(), xlim = c(0, 10)), class = "rnapeaks_error_invalid_arg")
  good <- data.frame(start = 0, end = 900, strand = "+", y = 0)
  expect_error(build_intron_arrows(good, xlim = c(0, 10, 20)), class = "rnapeaks_error_invalid_arg")
  expect_error(build_intron_arrows(good, xlim = c(10, 0)), class = "rnapeaks_error_invalid_arg")
  expect_error(build_intron_arrows(good, xlim = c(0, 900), min_intron_frac = 1.2),
               class = "rnapeaks_error_invalid_arg")
  expect_error(build_intron_arrows(good, xlim = c(0, 900), shaft_frac = -0.1),
               class = "rnapeaks_error_invalid_arg")
  expect_error(build_intron_arrows(good, xlim = c(0, 900), arrows_per_view = -1),
               class = "rnapeaks_error_invalid_arg")
})
