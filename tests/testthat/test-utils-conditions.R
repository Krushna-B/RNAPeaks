# Tests for R/utils-conditions.R
#
# These assert the *documented contract* of each helper, not the current
# implementation's incidental output. A failing test here means the helper
# does not do what its doc comment promises.
#
# Error-class contract (from .make_aborter / .abort_check):
#   - every RNAPeaks abort carries the umbrella class "rnapeaks_error"
#   - validators (check_*) abort with "rnapeaks_error_invalid_arg"

# --- classed aborters: correct class hierarchy ----------------------------

test_that("abort_* helpers carry their specific class plus rnapeaks_error", {
  expect_error(abort_invalid_bed("bad"), class = "rnapeaks_error_invalid_bed")
  expect_error(abort_invalid_bed("bad"), class = "rnapeaks_error")

  expect_error(abort_invalid_arg("bad"), class = "rnapeaks_error_invalid_arg")
  expect_error(abort_invalid_arg("bad"), class = "rnapeaks_error")

  expect_error(abort_not_found("bad"), class = "rnapeaks_error_not_found")
  expect_error(abort_not_found("bad"), class = "rnapeaks_error")

  expect_error(abort_invalid_gtf("bad"), class = "rnapeaks_error_invalid_gtf")
  expect_error(abort_invalid_gtf("bad"), class = "rnapeaks_error")
})

test_that("abort_species_unknown is a specialised invalid_gtf error", {
  expect_error(abort_species_unknown(), class = "rnapeaks_error_species_unknown")
  expect_error(abort_species_unknown(), class = "rnapeaks_error_invalid_gtf")
  expect_error(abort_species_unknown(), class = "rnapeaks_error")
})

# --- check_scalar_number --------------------------------------------------

test_that("check_scalar_number accepts a single in-range number and returns it invisibly", {
  expect_equal(check_scalar_number(5, "a"), 5)
  expect_equal(check_scalar_number(0.5, "a", min = 0, max = 1), 0.5)
  expect_invisible(check_scalar_number(5, "a"))
})

test_that("check_scalar_number accepts the inclusive boundaries", {
  expect_equal(check_scalar_number(0, "a", min = 0, max = 10), 0)
  expect_equal(check_scalar_number(10, "a", min = 0, max = 10), 10)
})

test_that("check_scalar_number rejects non-numeric, wrong length, NA, and out-of-range", {
  expect_error(check_scalar_number("5", "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_scalar_number(c(1, 2), "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_scalar_number(numeric(0), "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_scalar_number(NA_real_, "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_scalar_number(11, "a", min = 0, max = 10), class = "rnapeaks_error_invalid_arg")
  expect_error(check_scalar_number(-1, "a", min = 0, max = 10), class = "rnapeaks_error_invalid_arg")
})

# --- check_scalar_int -----------------------------------------------------

test_that("check_scalar_int accepts whole numbers (integer or double-valued)", {
  expect_equal(check_scalar_int(3L, "a"), 3L)
  expect_equal(check_scalar_int(3, "a"), 3)
  expect_equal(check_scalar_int(0, "a", min = 0, max = 5), 0)
})

test_that("check_scalar_int rejects fractional values", {
  expect_error(check_scalar_int(1.5, "a"), class = "rnapeaks_error_invalid_arg")
})

test_that("check_scalar_int rejects non-numeric, wrong length, NA, and out-of-range", {
  expect_error(check_scalar_int("3", "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_scalar_int(c(1L, 2L), "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_scalar_int(NA_integer_, "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_scalar_int(6, "a", min = 0, max = 5), class = "rnapeaks_error_invalid_arg")
})

# --- check_unit_interval --------------------------------------------------

test_that("check_unit_interval accepts [0, 1] and rejects outside", {
  expect_equal(check_unit_interval(0, "a"), 0)
  expect_equal(check_unit_interval(1, "a"), 1)
  expect_equal(check_unit_interval(0.42, "a"), 0.42)
  expect_error(check_unit_interval(1.01, "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_unit_interval(-0.01, "a"), class = "rnapeaks_error_invalid_arg")
})

# --- check_scalar_number_or_null ------------------------------------------

test_that("check_scalar_number_or_null lets NULL through", {
  expect_null(check_scalar_number_or_null(NULL, "a"))
})

test_that("check_scalar_number_or_null validates non-NULL like check_scalar_number", {
  expect_equal(check_scalar_number_or_null(2, "a", min = 0, max = 3), 2)
  expect_error(check_scalar_number_or_null(9, "a", min = 0, max = 3), class = "rnapeaks_error_invalid_arg")
  expect_error(check_scalar_number_or_null(NA_real_, "a"), class = "rnapeaks_error_invalid_arg")
})

# --- check_flag -----------------------------------------------------------

test_that("check_flag accepts a single TRUE/FALSE", {
  expect_true(check_flag(TRUE, "a"))
  expect_false(check_flag(FALSE, "a"))
})

test_that("check_flag rejects non-logical, wrong length, and NA", {
  expect_error(check_flag(1, "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_flag("TRUE", "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_flag(c(TRUE, FALSE), "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_flag(NA, "a"), class = "rnapeaks_error_invalid_arg")
})

# --- check_string ---------------------------------------------------------

test_that("check_string accepts a single non-NA string", {
  expect_equal(check_string("hi", "a"), "hi")
})

test_that("check_string enforces choices when supplied", {
  expect_equal(check_string("Count", "a", choices = c("Count", "Score")), "Count")
  expect_error(check_string("nope", "a", choices = c("Count", "Score")),
               class = "rnapeaks_error_invalid_arg")
})

test_that("check_string rejects non-character, wrong length, and NA", {
  expect_error(check_string(1, "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_string(c("a", "b"), "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_string(NA_character_, "a"), class = "rnapeaks_error_invalid_arg")
})

# --- check_color ----------------------------------------------------------

test_that("check_color accepts named colors and hex strings", {
  expect_equal(check_color("red", "a"), "red")
  expect_equal(check_color("#FF0000", "a"), "#FF0000")
})

test_that("check_color accepts NA only when allow_na = TRUE", {
  expect_true(is.na(check_color(NA, "a", allow_na = TRUE)))
  expect_error(check_color(NA, "a", allow_na = FALSE), class = "rnapeaks_error_invalid_arg")
})

test_that("check_color rejects unknown colors and non-scalar input", {
  expect_error(check_color("notacolor", "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_color(c("red", "blue"), "a"), class = "rnapeaks_error_invalid_arg")
})

# --- check_range_or_null --------------------------------------------------

test_that("check_range_or_null accepts NULL and an ascending length-2 numeric", {
  expect_null(check_range_or_null(NULL, "a"))
  expect_equal(check_range_or_null(c(0, 10), "a"), c(0, 10))
})

test_that("check_range_or_null rejects wrong length, NA, and min >= max", {
  expect_error(check_range_or_null(5, "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_range_or_null(c(1, 2, 3), "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_range_or_null(c(NA, 2), "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_range_or_null(c(10, 0), "a"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_range_or_null(c(5, 5), "a"), class = "rnapeaks_error_invalid_arg")
})

# --- %||% -----------------------------------------------------------------

test_that("%||% returns the left side unless it is NULL", {
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% 2, 2)
  expect_equal(NA %||% 2, NA)   # NA is not NULL, so left side wins
})

# --- normalize_str --------------------------------------------------------

test_that("normalize_str trims a single string and passes NULL through", {
  expect_equal(normalize_str("  hi  "), "hi")
  expect_null(normalize_str(NULL))
})

test_that("normalize_str leaves non-scalar-string input unchanged", {
  expect_equal(normalize_str(5), 5)
  expect_equal(normalize_str(c("a ", "b ")), c("a ", "b "))
  expect_true(is.na(normalize_str(NA_character_)))
})

# --- normalize_chr --------------------------------------------------------

test_that("normalize_chr strips a leading CHR (any case), trims, and uppercases", {
  expect_equal(normalize_chr("chr19"), "19")
  expect_equal(normalize_chr("Chr19"), "19")
  expect_equal(normalize_chr(" CHR19 "), "19")
  expect_equal(normalize_chr("chrX"), "X")
  expect_equal(normalize_chr("19"), "19")
  expect_equal(normalize_chr(19), "19")
})

test_that("normalize_chr passes NULL through", {
  expect_null(normalize_chr(NULL))
})

# --- normalize_coord ------------------------------------------------------

test_that("normalize_coord returns numeric unchanged and passes NULL through", {
  expect_equal(normalize_coord(43190342), 43190342)
  expect_null(normalize_coord(NULL))
})

test_that("normalize_coord parses comma-grouped and padded numeric strings", {
  expect_equal(normalize_coord("43,190,342"), 43190342)
  expect_equal(normalize_coord("  1200 "), 1200)
  expect_equal(normalize_coord(c("1", "2,000")), c(1, 2000))
})

test_that("normalize_coord aborts on unparseable strings", {
  expect_error(normalize_coord("abc"), class = "rnapeaks_error_invalid_arg")
  expect_error(normalize_coord(c("1", "oops")), class = "rnapeaks_error_invalid_arg")
})

# --- check_data_frame -----------------------------------------------------

test_that("check_data_frame accepts a non-empty frame containing required columns", {
  df <- data.frame(chr = "1", start = 1L, end = 2L)
  expect_equal(check_data_frame(df, "df", required = c("chr", "start")), df)
})

test_that("check_data_frame rejects non-frames, empty frames, and missing columns", {
  expect_error(check_data_frame(list(a = 1), "df"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_data_frame(data.frame(chr = character(0)), "df"),
               class = "rnapeaks_error_invalid_arg")
  expect_error(
    check_data_frame(data.frame(chr = "1"), "df", required = c("chr", "start")),
    class = "rnapeaks_error_invalid_arg"
  )
})
