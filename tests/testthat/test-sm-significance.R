# Tests for the per-position significance test in R/sm-event_significance.R
#
# Both stat_test modes run the SAME one-sided hypergeometric (Fisher) test; they
# differ only in the Control counts fed to it. So the p-values are exact and
# hand-derivable rather than stochastic. Expected values below are worked out
# from the hypergeometric mass function directly, not read back from phyper():
#
#   p = P(X >= group_hits),  X ~ Hypergeometric(total_hits, total_misses, n_group)
#     = sum_{x >= group_hits} C(total_hits, x) C(total_misses, n_group - x)
#                             / C(total_hits + total_misses, n_group)
#
# Worked cases (n_group = n_control = 5, so C(10, 5) = 252):
#   group 3 / control 1: total_hits 4, misses 6 -> (60 + 6)/252 = 66/252 = 11/42
#   group 5 / control 0: total_hits 5, misses 5 -> 1/252
#   group 0 / anything : P(X >= 0) = 1

test_that(".fisher_per_position matches the hand-derived hypergeometric tail", {
  p <- .fisher_per_position(
    group_hits   = c(3, 5, 0),
    n_group      = 5,
    control_hits = c(1, 0, 2),
    n_control    = 5
  )
  expect_equal(p, c(11 / 42, 1 / 252, 1))
})

test_that(".fisher_per_position is monotone: more group hits -> smaller p", {
  # Hold the Control cell fixed and walk group_hits up; enrichment can only grow.
  p <- .fisher_per_position(group_hits = 0:5, n_group = 5,
                            control_hits = rep(1, 6), n_control = 5)
  expect_true(all(diff(p) <= 0))
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(p[1], 1)                 # zero group hits is never significant
})

test_that(".fisher_per_position is vectorised over positions", {
  p <- .fisher_per_position(group_hits = c(1, 2), n_group = 5,
                            control_hits = c(0, 0), n_control = 5)
  expect_length(p, 2)
})

# --- .bootstrap_control_cell ----------------------------------------------
# Rebuilds a whole-count Control cell from the bootstrap's per-position hit
# frequency: counts = round(mean * sample_size), n_control = sample_size.

test_that(".bootstrap_control_cell rounds mean * sample_size to whole counts", {
  stats <- list(mean_per_position = c(0.2, 0.5, 0.44), sample_size = 10)
  cell  <- .bootstrap_control_cell(stats, control_counts = c(9, 9, 9),
                                   n_control = 99)
  expect_equal(cell$control_counts, c(2, 5, 4))   # round(c(2, 5, 4.4))
  expect_identical(cell$n_control, 10L)
})

test_that(".bootstrap_control_cell falls back to the raw pool when degenerate", {
  raw_counts <- c(3, 1)
  # NULL stats and a sub-1 sample size both mean "no usable bootstrap": warn and
  # hand back exactly the raw Control cell that was passed in.
  expect_warning(
    null_cell <- .bootstrap_control_cell(NULL, raw_counts, n_control = 5),
    regexp = "raw Control"
  )
  expect_equal(null_cell$control_counts, raw_counts)
  expect_equal(null_cell$n_control, 5)

  expect_warning(
    tiny_cell <- .bootstrap_control_cell(
      list(mean_per_position = c(0.5, 0.5), sample_size = 0),
      raw_counts, n_control = 5
    ),
    regexp = "raw Control"
  )
  expect_equal(tiny_cell$control_counts, raw_counts)
})

# --- test_per_position ----------------------------------------------------

test_that("fisher-all builds the p-value / FDR / significance frame", {
  # Raw p per position (from the worked cases above): 1/252, 11/42, 1.
  opts <- splicing_options(stat_test = "fisher-all", use_fdr = FALSE,
                           fdr_threshold = 0.05)
  res <- test_per_position(
    group_counts   = c(5, 3, 0), n_group = 5,
    control_counts = c(0, 1, 2), n_control = 5,
    opts = opts
  )
  expect_equal(res$position, 1:3)
  expect_equal(res$pvalue, c(1 / 252, 11 / 42, 1))
  expect_equal(res$pvalue_adj, res$pvalue)          # no FDR -> adj == raw
  expect_equal(res$significant, c(TRUE, FALSE, FALSE))
})

test_that("use_fdr applies Benjamini-Hochberg to the raw p-values", {
  # BH on sorted (1/252, 11/42, 1): adj = (p * n / rank), enforced monotone.
  #   rank 1: (1/252) * 3/1 = 3/252 = 1/84
  #   rank 2: (11/42) * 3/2 = 33/84 = 11/28
  #   rank 3: 1
  opts <- splicing_options(stat_test = "fisher-all", use_fdr = TRUE,
                           fdr_threshold = 0.05)
  res <- test_per_position(
    group_counts   = c(5, 3, 0), n_group = 5,
    control_counts = c(0, 1, 2), n_control = 5,
    opts = opts
  )
  expect_equal(res$pvalue_adj, c(1 / 84, 11 / 28, 1))
  expect_equal(res$significant, c(TRUE, FALSE, FALSE))
})

test_that("fisher-bootstrap substitutes the resampled Control cell", {
  # Bootstrap says the Control rate is 0 everywhere, so the substituted cell has
  # 0 hits regardless of the raw Control counts. group 5 / control 0 -> 1/252,
  # which fisher-all (control 3) would NOT produce -- proving the substitution.
  opts  <- splicing_options(stat_test = "fisher-bootstrap", use_fdr = FALSE)
  stats <- list(mean_per_position = 0, sample_size = 5)
  res <- test_per_position(
    group_counts   = 5, n_group = 5,
    control_counts = 3, n_control = 5,
    control_stats  = stats, opts = opts
  )
  expect_equal(res$pvalue, 1 / 252)
})

test_that("fisher-bootstrap warns and uses raw Control when the bootstrap is unusable", {
  opts <- splicing_options(stat_test = "fisher-bootstrap", use_fdr = FALSE)
  expect_warning(
    res <- test_per_position(
      group_counts   = 5, n_group = 5,
      control_counts = 0, n_control = 5,
      control_stats  = NULL, opts = opts
    ),
    regexp = "raw Control"
  )
  expect_equal(res$pvalue, 1 / 252)   # raw control 0 -> same tail as above
})
