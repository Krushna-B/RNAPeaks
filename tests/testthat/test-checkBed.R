test_that("checkBed validates proper BED format", {
  # Valid BED data
  valid_bed <- data.frame(
    chr = c("chr1", "chr1"),
    start = c(100, 200),
    end = c(150, 250),
    name = c("peak1", "peak2"),
    score = c(100, 200),
    strand = c("+", "+"),
    stringsAsFactors = FALSE
  )

  result <- checkBed(valid_bed)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true("chr" %in% colnames(result))
  expect_true("start" %in% colnames(result))
  expect_true("end" %in% colnames(result))
  # chr prefix should be removed
  expect_equal(result$chr[1], "1")
})

test_that("checkBed removes chr prefix", {
  bed <- data.frame(
    chr = c("chr17", "chrX"),
    start = c(100, 200),
    end = c(150, 250),
    name = c("peak1", "peak2"),
    score = c(100, 200),
    strand = c("+", "-"),
    stringsAsFactors = FALSE
  )

  result <- checkBed(bed)

  expect_equal(result$chr[1], "17")
  expect_equal(result$chr[2], "x")
})

test_that("checkBed errors on invalid strand", {
  invalid_bed <- data.frame(
    chr = c("chr1"),
    start = c(100),
    end = c(150),
    name = c("peak1"),
    score = c(100),
    strand = c("invalid"),
    stringsAsFactors = FALSE
  )

  expect_error(checkBed(invalid_bed), "Check Bed File!")
})

test_that("checkBed errors on end < start", {
  invalid_bed <- data.frame(
    chr = c("chr1"),
    start = c(200),
    end = c(100),
    name = c("peak1"),
    score = c(100),
    strand = c("+"),
    stringsAsFactors = FALSE
  )

  expect_error(checkBed(invalid_bed), "Check Bed File!")
})

test_that("checkBed requires at least 3 columns", {
  invalid_bed <- data.frame(
    chr = c("chr1"),
    start = c(100),
    stringsAsFactors = FALSE
  )

  expect_error(checkBed(invalid_bed))
})
