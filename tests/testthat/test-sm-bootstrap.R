# Tests for bootstrap_control() in R/sm-event_bootstrap.R
#
# The bootstrap is stochastic, so instead of re-deriving the random draws we
# anchor on properties that hold for ANY correct implementation:
#   * a position hit by ALL controls has resample rate 1 with zero SD, and a
#     position hit by NONE has rate 0 -- true regardless of the RNG;
#   * the bootstrap mean must CONVERGE to the empirical control rate (the whole
#     point of the procedure) -- checked to 1/3 for a 1-of-3 position;
#   * exact analytic branches (sample_size < 1, no controls) and reproducibility.
#
# Fixture control matrix (n_control = 3, n_positions = 4):
#   pos 1: hit by all 3 controls  -> rate 1
#   pos 2: hit by none            -> rate 0
#   pos 3: hit by control 1 only  -> rate 1/3
#   pos 4: hit by none            -> rate 0
ctrl_hits <- list(event_id = c(1L, 2L, 3L, 1L), col_idx = c(1L, 1L, 1L, 3L))

test_that("fully-covered and empty positions are exact regardless of the seed", {
  opts <- splicing_options(control_multiplier = 2, control_iterations = 50)
  set.seed(1)
  res <- suppressMessages(
    bootstrap_control(ctrl_hits, n_control = 3, n_positions = 4,
                      n_pos = 2, n_neg = 0, opts = opts)
  )
  expect_equal(res$mean_per_position[c(1, 2, 4)], c(1, 0, 0))
  expect_equal(res$sd_per_position[c(1, 2, 4)], c(0, 0, 0))
  # a partially-covered position must carry some resampling variance
  expect_gt(res$sd_per_position[3], 0)
  expect_true(res$mean_per_position[3] >= 0 && res$mean_per_position[3] <= 1)
})

test_that("the bootstrap mean converges to the empirical control rate", {
  opts <- splicing_options(control_multiplier = 2, control_iterations = 3000)
  set.seed(42)
  res <- suppressMessages(
    bootstrap_control(ctrl_hits, 3, 4, n_pos = 2, n_neg = 0, opts = opts)
  )
  # position 3 is hit by 1 of 3 controls -> mean must approach 1/3
  expect_lt(abs(res$mean_per_position[3] - 1 / 3), 0.03)
  expect_equal(res$mean_per_position[1], 1)   # 3 of 3
  expect_equal(res$mean_per_position[2], 0)   # 0 of 3
})

test_that("bootstrap_control is reproducible under a fixed seed", {
  opts <- splicing_options(control_multiplier = 2, control_iterations = 40)
  set.seed(7); r1 <- suppressMessages(bootstrap_control(ctrl_hits, 3, 4, 2, 0, opts))
  set.seed(7); r2 <- suppressMessages(bootstrap_control(ctrl_hits, 3, 4, 2, 0, opts))
  expect_equal(r1, r2)
})

# --- analytic branches ----------------------------------------------------

test_that("sample_size < 1 returns the raw control mean with zero SD", {
  opts <- splicing_options(control_multiplier = 2, control_iterations = 20)
  res <- suppressWarnings(
    bootstrap_control(ctrl_hits, 3, 4, n_pos = 0, n_neg = 0, opts = opts)
  )
  expect_equal(res$mean_per_position, c(1, 0, 1 / 3, 0))   # counts / n_control
  expect_equal(res$sd_per_position, rep(0, 4))
})

test_that("no controls warns and returns a flat zero mean and SD", {
  opts <- splicing_options()
  expect_warning(bootstrap_control(NULL, 0, 4, 2, 0, opts), regexp = "No Control")
  res <- suppressWarnings(bootstrap_control(NULL, 0, 4, 2, 0, opts))
  expect_equal(res$mean_per_position, rep(0, 4))
  expect_equal(res$sd_per_position, rep(0, 4))
})

test_that("a non-null sample size but empty control hits gives zeros", {
  opts <- splicing_options(control_multiplier = 2, control_iterations = 20)
  res <- suppressMessages(bootstrap_control(NULL, 3, 4, 2, 0, opts))
  expect_equal(res$mean_per_position, rep(0, 4))
})

test_that("bootstrap_control validates n_pos / n_neg", {
  opts <- splicing_options()
  expect_error(bootstrap_control(ctrl_hits, 3, 4, n_pos = -1, n_neg = 0, opts),
               class = "rnapeaks_error_invalid_arg")
  expect_error(bootstrap_control(ctrl_hits, 3, 4, n_pos = 0, n_neg = -1, opts),
               class = "rnapeaks_error_invalid_arg")
})
