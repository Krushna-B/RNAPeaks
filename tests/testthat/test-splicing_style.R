# Tests for splicing_style() in R/params-splicing_style.R
#
# Contract: validated named list of visual settings. group_colors must supply
# non-NA colors for Positive / Negative / Control; other args delegate to the
# check_* validators.

test_that("defaults produce a valid list with documented values", {
  st <- splicing_style()
  expect_type(st, "list")
  expect_equal(st$legend_position, "bottom")
  expect_equal(st$ylab, "Frequency")
})

test_that("the returned list contains exactly the documented parameters", {
  # Surfaces any internal temp variables leaking out of as.list(environment()).
  expect_setequal(names(splicing_style()), names(formals(splicing_style)))
})

test_that("group_colors must be a named list covering all three groups with valid colors", {
  expect_error(splicing_style(group_colors = c("blue", "red", "black")),
               class = "rnapeaks_error_invalid_arg")
  expect_error(
    splicing_style(group_colors = list(Positive = "blue", Negative = "red")),
    class = "rnapeaks_error_invalid_arg"
  )
  expect_error(
    splicing_style(group_colors = list(Positive = "notacolor", Negative = "red", Control = "black")),
    class = "rnapeaks_error_invalid_arg"
  )
  expect_error(
    splicing_style(group_colors = list(Positive = NA, Negative = "red", Control = "black")),
    class = "rnapeaks_error_invalid_arg"
  )
})

test_that("line/ribbon numeric ranges are enforced", {
  expect_error(splicing_style(line_width = -1), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_style(line_alpha = 1.5), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_style(ribbon_alpha = -0.1), class = "rnapeaks_error_invalid_arg")
})

test_that("flags, colors and legend choices are validated", {
  expect_error(splicing_style(show_significance = NA), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_style(title_color = "notacolor"), class = "rnapeaks_error_invalid_arg")
  expect_error(splicing_style(legend_position = "middle"), class = "rnapeaks_error_invalid_arg")
})

test_that("isoform label nudges accept negative values but sizes do not", {
  expect_equal(splicing_style(isoform_label_nudge_x = -5)$isoform_label_nudge_x, -5)
  expect_error(splicing_style(isoform_label_size = -1), class = "rnapeaks_error_invalid_arg")
})

# --- exhaustive: every argument's validation is wired up ------------------

test_that("every argument rejects a known-bad value with an invalid_arg error", {
  bad_values <- list(
    group_colors          = c("blue", "red", "black"),  # not a named list
    line_width            = -1,
    line_alpha            = 1.5,
    ribbon_alpha          = 1.5,
    show_significance     = NA,
    title_size            = -1,
    title_color           = "notacolor",
    axis_text_size        = -1,
    ylab                  = 1,          # must be a string
    boundary_col          = "notacolor",
    exon_col              = "notacolor",
    isoform_label_size    = -1,
    isoform_label_nudge_x = "x",        # free numeric -> reject non-numeric
    isoform_label_nudge_y = "x",
    legend_position       = "middle"
  )
  expect_setequal(names(bad_values), names(formals(splicing_style)))
  for (arg in names(bad_values)) {
    expect_error(
      do.call(splicing_style, stats::setNames(list(bad_values[[arg]]), arg)),
      class = "rnapeaks_error_invalid_arg",
      info  = arg
    )
  }
})
