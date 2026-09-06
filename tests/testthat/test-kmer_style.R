# Tests for kmer_style() in R/params-kmer_style.R
#
# Contract: validated named list of visual settings for the kmer_enrichment()
# scatter / rank plots. Colors reject NA / invalid; sizes must be non-negative;
# point_alpha is a unit interval; ref_line_type is a fixed choice;
# label_check_overlap is a flag.

test_that("defaults produce a valid list with documented values", {
  st <- kmer_style()
  expect_type(st, "list")
  expect_equal(st$point_color, "#3B6FB6")
  expect_equal(st$ref_line_type, "dashed")
  expect_true(st$label_check_overlap)
})

test_that("the returned list contains exactly the documented parameters", {
  expect_setequal(names(kmer_style()), names(formals(kmer_style)))
})

test_that("ref_line_type is limited to the ggplot2 line types", {
  expect_equal(kmer_style(ref_line_type = "dotted")$ref_line_type, "dotted")
  expect_error(kmer_style(ref_line_type = "wavy"),
               class = "rnapeaks_error_invalid_arg")
})

test_that("colors reject NA and invalid values", {
  expect_error(kmer_style(point_color = "notacolor"),
               class = "rnapeaks_error_invalid_arg")
  expect_error(kmer_style(ref_line_color = NA),
               class = "rnapeaks_error_invalid_arg")
  expect_error(kmer_style(title_color = "notacolor"),
               class = "rnapeaks_error_invalid_arg")
})

# --- exhaustive: every argument's validation is wired up ------------------

test_that("every argument rejects a known-bad value with an invalid_arg error", {
  bad_values <- list(
    point_color         = "notacolor",
    point_size          = -1,
    point_alpha         = 1.5,
    label_size          = -1,
    label_check_overlap = "yes",        # must be a flag
    ref_line_color      = "notacolor",
    ref_line_type       = "wavy",
    title_size          = -1,
    title_color         = "notacolor",
    axis_text_size      = -1
  )
  expect_setequal(names(bad_values), names(formals(kmer_style)))
  for (arg in names(bad_values)) {
    expect_error(
      do.call(kmer_style, stats::setNames(list(bad_values[[arg]]), arg)),
      class = "rnapeaks_error_invalid_arg",
      info  = arg
    )
  }
})
