# Tests for check_bed() in R/bed-check.R
#
# Contract:
#   - Accepts a data frame OR a list of data frames, each with >= 6 columns.
#   - Columns are mapped by POSITION: 1=chr, 2=start, 3=end, 6=strand; those four
#     positions are renamed to canonical names, all other columns are preserved.
#   - Produces a single combined frame with a `target` column that groups rows
#     (by bed label and/or split_col value).
#   - Coordinates are coerced to numeric; rows with unparseable start/end are
#     dropped with a warning. Invariants enforced: start >= 0, end >= start,
#     strand in {"+","-"}.
#
# --- input shape ----------------------------------------------------------

test_that("a bare data frame is accepted and returns a single frame", {
  out <- check_bed(make_raw_bed(3))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3)
})

test_that("a list of frames is row-combined", {
  out <- check_bed(list(make_raw_bed(2), make_raw_bed(3)))
  expect_equal(nrow(out), 5)
})

test_that("non-frame / empty / non-list inputs abort as invalid_bed", {
  expect_error(check_bed("nope"), class = "rnapeaks_error_invalid_bed")
  expect_error(check_bed(list()), class = "rnapeaks_error_invalid_bed")
  expect_error(check_bed(list(1, 2)), class = "rnapeaks_error_invalid_bed")
  expect_error(check_bed(matrix(1:12, ncol = 6)), class = "rnapeaks_error_invalid_bed")
})

# --- column count ---------------------------------------------------------

test_that("fewer than 6 columns aborts", {
  expect_error(check_bed(make_raw_bed(2)[, 1:5]), class = "rnapeaks_error_invalid_bed")
})

test_that("beds with mismatched column counts cannot be combined", {
  wide <- cbind(make_raw_bed(2), extra = 1)
  expect_error(check_bed(list(make_raw_bed(2), wide)),
               class = "rnapeaks_error_invalid_bed")
})

# --- column naming & preservation -----------------------------------------

test_that("positions 1,2,3,6 are renamed to canonical regardless of input names", {
  bed <- make_raw_bed(2)
  colnames(bed) <- c("Chrom", "Begin", "Stop", "Prot", "Sig", "Str")
  out <- check_bed(bed)
  expect_equal(colnames(out)[c(1, 2, 3, 6)], c("chr", "start", "end", "strand"))
})

test_that("extra columns beyond position 6 are preserved by name", {
  bed <- cbind(make_raw_bed(2), extra = c("a", "b"))
  out <- check_bed(bed)
  expect_true("extra" %in% colnames(out))
})

test_that("beds that disagree on non-canonical column names abort", {
  b1 <- make_raw_bed(2)
  b2 <- make_raw_bed(2)
  colnames(b2)[5] <- "foo"          # position 5 is not canonicalised
  expect_error(check_bed(list(b1, b2)), class = "rnapeaks_error_invalid_bed")
})

# --- split_col validation -------------------------------------------------

test_that("split_col must be a positive whole number not pointing at a canonical column", {
  bed <- make_raw_bed(2)
  expect_error(check_bed(bed, split_col = 1), class = "rnapeaks_error_invalid_arg")
  expect_error(check_bed(bed, split_col = 6), class = "rnapeaks_error_invalid_arg")
  expect_error(check_bed(bed, split_col = 0), class = "rnapeaks_error_invalid_arg")
  expect_error(check_bed(bed, split_col = 1.5), class = "rnapeaks_error_invalid_arg")
  expect_error(check_bed(bed, split_col = "4"), class = "rnapeaks_error_invalid_arg")
  expect_error(check_bed(bed, split_col = NA), class = "rnapeaks_error_invalid_arg")
})

test_that("split_col exceeding the available columns aborts", {
  expect_error(check_bed(make_raw_bed(2), split_col = 7),
               class = "rnapeaks_error_invalid_arg")
})

# --- target column semantics ----------------------------------------------

test_that("single unnamed bed, no split_col -> target is the synthesised label", {
  out <- check_bed(make_raw_bed(2))
  expect_equal(unique(out$target), "bed1")
})

test_that("single named bed, no split_col -> target is the list name", {
  out <- check_bed(list(MyBed = make_raw_bed(2)))
  expect_equal(unique(out$target), "MyBed")
})

test_that("split_col draws target from that column's values", {
  out <- check_bed(make_raw_bed(3, protein = "SRSF1"), split_col = 4)
  expect_equal(unique(out$target), "SRSF1")
})

test_that("named multi-bed with split_col combines label and value", {
  out <- check_bed(
    list(K562 = make_raw_bed(2, protein = "SRSF1"),
         HepG2 = make_raw_bed(2, protein = "U2AF2")),
    split_col = 4
  )
  expect_setequal(unique(out$target), c("K562_SRSF1", "HepG2_U2AF2"))
})

test_that("unnamed multi-bed with split_col synthesises per-bed labels", {
  out <- check_bed(
    list(make_raw_bed(1, protein = "A"), make_raw_bed(1, protein = "B")),
    split_col = 4
  )
  expect_setequal(unique(out$target), c("bed1_A", "bed2_B"))
})

test_that("named multi-bed, no split_col -> target is the label alone", {
  out <- check_bed(list(K562 = make_raw_bed(2), HepG2 = make_raw_bed(3)))
  expect_setequal(unique(out$target), c("K562", "HepG2"))
})

# --- coordinate normalisation & validation --------------------------------

test_that("chromosome names are normalised (chr prefix stripped, uppercased)", {
  out <- check_bed(make_raw_bed(2, chr = "chr7"))
  expect_equal(unique(out$chr), "7")
})

test_that("character/factor coordinates are coerced to numeric", {
  bed <- make_raw_bed(2)
  bed$start <- as.character(bed$start)
  bed$end <- factor(as.character(bed$end))
  out <- check_bed(bed)
  expect_type(out$start, "double")
  expect_type(out$end, "double")
})

test_that("rows with unparseable coordinates are dropped with a warning", {
  bed <- make_raw_bed(3)
  bed$start <- c("100", "not_a_number", "500")
  expect_warning(check_bed(bed), regexp = "[Dd]ropped")
  out <- suppressWarnings(check_bed(bed))
  expect_equal(nrow(out), 2)
})

test_that("a bed whose coordinates are all unparseable aborts", {
  bed <- make_raw_bed(2)
  bed$start <- c("x", "y")
  expect_error(suppressWarnings(check_bed(bed)), class = "rnapeaks_error_invalid_bed")
})

test_that("negative coordinates abort", {
  expect_error(check_bed(make_raw_bed(2, start = c(-5, 100), end = c(50, 150))),
               class = "rnapeaks_error_invalid_bed")
  expect_error(check_bed(make_raw_bed(1, start = 100, end = -5)),
               class = "rnapeaks_error_invalid_bed")
})

test_that("end < start aborts, but a zero-width peak (end == start) is allowed", {
  expect_error(check_bed(make_raw_bed(1, start = 200, end = 100)),
               class = "rnapeaks_error_invalid_bed")
  expect_silent(out <- check_bed(make_raw_bed(1, start = 100, end = 100)))
  expect_equal(nrow(out), 1)
})

test_that("start == 0 (0-based) is accepted", {
  out <- check_bed(make_raw_bed(1, start = 0, end = 50))
  expect_equal(nrow(out), 1)
})

# --- strand validation ----------------------------------------------------

test_that("strand must be + or -; anything else aborts", {
  expect_error(check_bed(make_raw_bed(2, strand = ".")),
               class = "rnapeaks_error_invalid_bed")
  expect_error(check_bed(make_raw_bed(2, strand = "*")),
               class = "rnapeaks_error_invalid_bed")
})

test_that("a mix of + and - strands is accepted", {
  bed <- make_raw_bed(2)
  bed$strand <- c("+", "-")
  out <- check_bed(bed)
  expect_setequal(out$strand, c("+", "-"))
})
