# Tests for splicing_options() in R/params-splicing_options.R
#
# Contract: validated named list of analysis options.

test_that("defaults produce a valid list with documented values", {
  opt <- splicing_options()
  expect_type(opt, "list")
  expect_equal(opt$psi_cutoff, c(-0.1, 0.1))
  expect_equal(opt$stat_test, "fisher-all")
  expect_true(opt$use_fdr)
})

test_that("the returned list contains exactly the documented parameters", {
  # Surfaces any internal temp variables leaking out of as.list(environment()).
  expect_setequal(names(splicing_options()), names(formals(splicing_options)))
})

# --- region geometry & simple ranges --------------------------------------

test_that("width_exon / width_intron must be positive integers", {
  expect_error(splicing_options(width_exon = 0), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(width_intron = 1.5), class = "rnapeaks_error_invalid_arg")
})

test_that("event_fdr and control_pval are unit-interval values", {
  expect_error(splicing_options(event_fdr = 1.2), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(control_pval = -0.1), class = "rnapeaks_error_invalid_arg")
})

# --- psi_cutoff invariants ------------------------------------------------

test_that("psi_cutoff must be a length-2 numeric within [-1, 1] and increasing", {
  expect_error(splicing_options(psi_cutoff = 0.1), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(psi_cutoff = c(-0.1, 0.1, 0.2)), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(psi_cutoff = c(-2, 0.1)), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(psi_cutoff = c(0.1, -0.1)), class = "rnapeaks_error_invalid_arg")  # neg >= pos
  expect_error(splicing_options(psi_cutoff = c(NA, 0.1)), class = "rnapeaks_error_invalid_arg")
})

test_that("psi_control_max must be strictly below min(abs(psi_cutoff))", {
  # default psi_cutoff = c(-0.1, 0.1) -> min(abs) = 0.1
  expect_error(splicing_options(psi_control_max = 0.1), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(psi_control_max = 0.2), class = "rnapeaks_error_invalid_arg")
  expect_silent(splicing_options(psi_control_max = 0.05))
})

# --- disable-via-NULL semantics -------------------------------------------

test_that("min_count = NULL disables the coverage filter (stored as 0)", {
  expect_equal(splicing_options(min_count = NULL)$min_count, 0L)
  expect_error(splicing_options(min_count = -1), class = "rnapeaks_error_invalid_arg")
})

test_that("moving_average = NULL disables smoothing (stored as 0)", {
  expect_equal(splicing_options(moving_average = NULL)$moving_average, 0L)
  expect_error(splicing_options(moving_average = -1), class = "rnapeaks_error_invalid_arg")
})

# --- groups ---------------------------------------------------------------

test_that("groups must be a non-empty, non-duplicated subset of the three groups", {
  expect_equal(splicing_options(groups = c("Positive", "Control"))$groups,
               c("Positive", "Control"))
  expect_error(splicing_options(groups = character(0)), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(groups = c("Positive", "Positive")), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(groups = c("Positive", "Foo")), class = "rnapeaks_error_invalid_arg")
})

# --- control bootstrap ----------------------------------------------------

test_that("control_multiplier must be strictly positive and iterations a positive integer", {
  expect_error(splicing_options(control_multiplier = 0), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(control_multiplier = -1), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(control_iterations = 0), class = "rnapeaks_error_invalid_arg")
})

# --- significance ---------------------------------------------------------

test_that("stat_test is limited to the fisher variants and flags/thresholds are validated", {
  expect_equal(splicing_options(stat_test = "fisher-all")$stat_test, "fisher-all")
  expect_equal(splicing_options(stat_test = "fisher-bootstrap")$stat_test, "fisher-bootstrap")
  expect_equal(splicing_options(stat_test = "fisher")$stat_test, "fisher-all")  # back-compat alias
  expect_error(splicing_options(stat_test = "binomial"), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(stat_test = "ttest"), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(use_fdr = "yes"), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(fdr_threshold = 2), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_options(verbose = NA), class = "rnapeaks_error_invalid_arg")
})

# --- exhaustive: every argument's validation is wired up ------------------
# min_count / moving_average accept NULL (disable) and so have no single-value
# rejection that fits the table; they are covered by dedicated tests above.

test_that("every argument rejects a known-bad value with an invalid_arg error", {
  bad_values <- list(
    width_exon         = 0,
    width_intron       = 0,
    event_fdr          = 1.2,
    control_pval       = -0.1,
    psi_cutoff         = c(0.1, -0.1),     # neg >= pos
    psi_control_max    = 0.2,              # >= min(abs(psi_cutoff))
    min_count          = -1,
    groups             = c("Positive", "Foo"),
    control_multiplier = 0,
    control_iterations = 0,
    moving_average     = -1,
    stat_test          = "ttest",
    use_fdr            = "yes",
    fdr_threshold      = 2,
    verbose            = NA
  )
  expect_setequal(names(bad_values), names(formals(splicing_options)))
  for (arg in names(bad_values)) {
    expect_error(
      do.call(splicing_options, stats::setNames(list(bad_values[[arg]]), arg)),
      class = "rnapeaks_error_invalid_arg",
      info  = arg
    )
  }
})
