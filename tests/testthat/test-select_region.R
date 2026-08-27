# Tests for select_region() in R/gene-select.R
#
# Contract: keep GTF rows for the LONGEST transcript of each gene whose body
# overlaps chr:start-end on the requested strand. Window overlap is by
# intersection (feature end >= start AND feature start <= end).

# --- validation -----------------------------------------------------------

test_that("select_region validates gtf shape and required columns", {
  expect_error(select_region(data.frame(), "1", 0, 100, "+"),
               class = "rnapeaks_error_invalid_arg")
  bad <- make_gtf()[, setdiff(names(make_gtf()), "width")]
  expect_error(select_region(bad, "1", 0, 2000, "+"),
               class = "rnapeaks_error_invalid_gtf")
})

test_that("select_region validates the window arguments", {
  gtf <- make_gtf()
  expect_error(select_region(gtf, 1, 0, 2000, "+"), class = "rnapeaks_error_invalid_arg")   # chr not string
  expect_error(select_region(gtf, "1", 0, 2000, "."), class = "rnapeaks_error_invalid_arg") # bad strand
  expect_error(select_region(gtf, "1", "x", 2000, "+"), class = "rnapeaks_error_invalid_arg")# start not numeric
  expect_error(select_region(gtf, "1", 2000, 0, "+"), class = "rnapeaks_error_invalid_arg")  # start > end
})

# --- per-gene longest-transcript selection --------------------------------

test_that("select_region keeps only the longest transcript of an overlapped gene", {
  out <- select_region(make_gtf(), "1", 0, 2000, "+")
  # window overlaps both SRSF1 transcripts; only the longest (ENST...1) is kept.
  expect_equal(unique(out$transcript_id), "ENST00000001")
  expect_true("transcript" %in% out$type)
})

test_that("chromosome prefix / case does not affect matching", {
  out <- select_region(make_gtf(), "chr1", 0, 2000, "+")
  expect_equal(unique(out$transcript_id), "ENST00000001")
})

test_that("strand is respected", {
  # 2000-2400 on + strand overlaps no gene; on - strand it hits MINUS.
  expect_error(select_region(make_gtf(), "1", 2000, 2400, "+"),
               class = "rnapeaks_error_not_found")
  out <- select_region(make_gtf(), "1", 2000, 2400, "-")
  expect_equal(unique(out$transcript_id), "ENST00000004")
})

test_that("a window spanning several genes returns one transcript per gene", {
  out <- select_region(make_gtf(), "1", 0, 6000, "+")
  expect_setequal(unique(out$transcript_id), c("ENST00000001", "ENST00000003"))
})

test_that("a window overlapping nothing aborts as not-found", {
  expect_error(select_region(make_gtf(), "2", 0, 100, "+"),
               class = "rnapeaks_error_not_found")
})
