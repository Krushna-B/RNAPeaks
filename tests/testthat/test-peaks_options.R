# Tests for peaks_options() in R/params-peaks_options.R
#
# Contract: return a validated named list of the documented peak-processing
# options. Invalid arguments abort with class "rnapeaks_error_invalid_arg".


test_that("defaults produce a valid list carrying the documented values", {
  opt <- peaks_options()
  expect_type(opt, "list")
  expect_equal(opt$split_by, 4)
  expect_equal(opt$order_by, "Count")
  expect_equal(opt$collapse, 0)
  expect_equal(opt$max_proteins, 100)
})

test_that("the returned list contains exactly the documented parameters", {
  expect_setequal(names(peaks_options()), names(formals(peaks_options)))
})

test_that("split_by is optional but otherwise a positive whole number", {
  expect_null(peaks_options(split_by = NULL)$split_by)
  expect_error(peaks_options(split_by = 0), class = "rnapeaks_error_invalid_arg")
  expect_error(peaks_options(split_by = 1.5), class = "rnapeaks_error_invalid_arg")
  expect_error(peaks_options(split_by = "4"), class = "rnapeaks_error_invalid_arg")
})

test_that("include / order_in accept NULL (off) or an NA-free character vector", {
  expect_null(peaks_options(include = NULL)$include)
  expect_equal(peaks_options(include = c("A", "B"))$include, c("A", "B"))
  expect_equal(peaks_options(order_in = c("A", "B"))$order_in, c("A", "B"))
})

test_that("include / order_in reject non-character or NA-containing vectors", {
  expect_error(peaks_options(include = 1L), class = "rnapeaks_error_invalid_arg")
  expect_error(peaks_options(include = c("A", NA)), class = "rnapeaks_error_invalid_arg")
  expect_error(peaks_options(order_in = c("A", NA)), class = "rnapeaks_error_invalid_arg")
})

test_that("order_by is limited to its two choices", {
  expect_equal(peaks_options(order_by = "Alphabetical")$order_by, "Alphabetical")
  expect_error(peaks_options(order_by = "Random"), class = "rnapeaks_error_invalid_arg")
})

test_that("collapse must be non-negative and max_proteins a positive integer", {
  expect_error(peaks_options(collapse = -1), class = "rnapeaks_error_invalid_arg")
  expect_error(peaks_options(max_proteins = 0), class = "rnapeaks_error_invalid_arg")
  expect_error(peaks_options(max_proteins = 1.5), class = "rnapeaks_error_invalid_arg")
})

# --- exhaustive: every argument's validation is wired up ------------------

test_that("every argument rejects a known-bad value with an invalid_arg error", {
  bad_values <- list(
    split_by     = 0,
    include      = 1L,             # must be character
    order_by     = "Random",
    order_in     = c("A", NA),     # NA not allowed
    collapse     = -1,
    max_proteins = 0
  )
  expect_setequal(names(bad_values), names(formals(peaks_options)))
  for (arg in names(bad_values)) {
    expect_error(
      do.call(peaks_options, stats::setNames(list(bad_values[[arg]]), arg)),
      class = "rnapeaks_error_invalid_arg",
      info  = arg
    )
  }
})
