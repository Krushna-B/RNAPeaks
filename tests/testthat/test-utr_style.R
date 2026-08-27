# Tests for utr_style() in R/params-utr_style.R
#
# Contract: validated named list of visual settings for plot_utr_binding().
# palette must be a non-empty NA-free character vector of valid colors;
# label_position is limited to "ends" / "center".

test_that("defaults produce a valid list with documented values", {
  st <- utr_style()
  expect_type(st, "list")
  expect_equal(st$label_position, "ends")
  expect_equal(st$single_track_color, "blue")
})

test_that("the returned list contains exactly the documented parameters", {
  # Surfaces any internal temp variables leaking out of as.list(environment()).
  expect_setequal(names(utr_style()), names(formals(utr_style)))
})

test_that("palette must be a non-empty vector of valid colors", {
  expect_equal(utr_style(palette = c("red", "blue"))$palette, c("red", "blue"))
  expect_error(utr_style(palette = character(0)), class = "rnapeaks_error_invalid_arg")
  expect_error(utr_style(palette = c("red", NA)), class = "rnapeaks_error_invalid_arg")
  expect_error(utr_style(palette = c("red", "notacolor")), class = "rnapeaks_error_invalid_arg")
})

test_that("label_position is limited to its two choices", {
  expect_equal(utr_style(label_position = "center")$label_position, "center")
  expect_error(utr_style(label_position = "top"), class = "rnapeaks_error_invalid_arg")
})

test_that("schematic colors and single_track_color reject NA / invalid colors", {
  expect_error(utr_style(utr_fill = NA), class = "rnapeaks_error_invalid_arg")
  expect_error(utr_style(cds_fill = "notacolor"), class = "rnapeaks_error_invalid_arg")
  expect_error(utr_style(single_track_color = NA), class = "rnapeaks_error_invalid_arg")
})

test_that("numeric sizes must be non-negative and legend choices are enforced", {
  expect_error(utr_style(line_width = -1), class = "rnapeaks_error_invalid_arg")
  expect_error(utr_style(schematic_height = -0.1), class = "rnapeaks_error_invalid_arg")
  expect_error(utr_style(line_alpha = 1.5), class = "rnapeaks_error_invalid_arg")
  expect_error(utr_style(legend_position = "middle"), class = "rnapeaks_error_invalid_arg")
})

# --- exhaustive: every argument's validation is wired up ------------------

test_that("every argument rejects a known-bad value with an invalid_arg error", {
  bad_values <- list(
    line_width         = -1,
    line_alpha         = 1.5,
    utr_fill           = "notacolor",
    cds_fill           = "notacolor",
    schematic_height   = -1,
    label_size         = -1,
    label_position     = "top",
    pct_label_size     = -1,
    pct_label_color    = "notacolor",
    palette            = character(0),
    single_track_color = "notacolor",
    title_size         = -1,
    title_color        = "notacolor",
    axis_text_size     = -1,
    ylab               = 1,          # must be a string
    legend_position    = "middle"
  )
  expect_setequal(names(bad_values), names(formals(utr_style)))
  for (arg in names(bad_values)) {
    expect_error(
      do.call(utr_style, stats::setNames(list(bad_values[[arg]]), arg)),
      class = "rnapeaks_error_invalid_arg",
      info  = arg
    )
  }
})
