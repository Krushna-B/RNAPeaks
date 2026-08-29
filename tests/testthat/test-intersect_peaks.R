# Tests for intersect_peaks() and its helpers in R/intersect_peaks.R
#
# intersect_peaks is an R reimplementation of `bedtools intersect -wa -wb`, so
# the overlap results are exact and hand-derivable. Inputs are supplied as data
# frames (BED layout, 0-based start) to keep the tests file-free; the helpers
# that parse coordinates / strand / GRanges are checked directly.

# BED6 row builder (start is 0-based, as in a real BED file).
bed6 <- function(chr, start, end, strand, name = "x", score = ".") {
  data.frame(chr = chr, start = start, end = end, name = name,
             score = score, strand = strand, stringsAsFactors = FALSE)
}

# --- .coerce_coord --------------------------------------------------------

test_that(".coerce_coord parses non-negative integers and rejects the rest", {
  expect_equal(.coerce_coord(c("10", "20"), "peaks", "start"), c(10L, 20L))
  expect_error(.coerce_coord(c("10", "x"), "peaks", "start"),
               class = "rnapeaks_error_invalid_bed")
  expect_error(.coerce_coord(c("-1", "5"), "peaks", "start"),
               class = "rnapeaks_error_invalid_bed")
})

# --- .coerce_strand -------------------------------------------------------

test_that(".coerce_strand maps '.' and NA to '*' and rejects unknown symbols", {
  expect_equal(.coerce_strand(c("+", "-", ".", NA), "peaks"),
               c("+", "-", "*", "*"))
  expect_error(.coerce_strand("?", "peaks"),
               class = "rnapeaks_error_invalid_bed")
})

# --- .bed_to_gr -----------------------------------------------------------

test_that(".bed_to_gr shifts the 0-based start to a 1-based GRanges", {
  gr <- .bed_to_gr(bed6("chr1", "100", "200", "+"), "peaks")
  expect_equal(GenomicRanges::start(gr), 101)   # 0-based 100 -> 1-based 101
  expect_equal(GenomicRanges::end(gr), 200)
  expect_equal(as.character(GenomicRanges::strand(gr)), "+")
})

test_that(".bed_to_gr rejects end < start", {
  expect_error(.bed_to_gr(bed6("chr1", "200", "100", "+"), "peaks"),
               class = "rnapeaks_error_invalid_bed")
})

# --- .read_bed_df ---------------------------------------------------------

test_that(".read_bed_df coerces a data frame to character and keeps names", {
  df <- .read_bed_df(bed6("chr1", 100, 200, "+"), "peaks")
  expect_type(df$start, "character")
  expect_equal(names(df), c("chr", "start", "end", "name", "score", "strand"))
})

test_that(".read_bed_df rejects an empty frame and a missing file", {
  expect_error(.read_bed_df(bed6("chr1", 100, 200, "+")[0, ], "peaks"),
               class = "rnapeaks_error_invalid_bed")
  expect_error(.read_bed_df("/no/such/file.bed", "peaks"),
               class = "rnapeaks_error_not_found")
})

# --- .normalize_peaks / .normalize_transcripts ----------------------------

test_that(".normalize_peaks requires at least six columns", {
  expect_error(.normalize_peaks(data.frame(a = 1, b = 2)),
               class = "rnapeaks_error_invalid_bed")
  expect_error(.normalize_peaks(42), class = "rnapeaks_error_invalid_arg")
})

test_that(".normalize_transcripts converts a GTF frame to 0-based BED6", {
  gtf <- data.frame(
    seqnames = "chr1", start = c(101, 101), end = c(200, 150),
    strand = "+", type = c("transcript", "exon"),
    transcript_id = "ENST1", stringsAsFactors = FALSE
  )
  bed <- .normalize_transcripts(gtf)
  expect_equal(nrow(bed), 1L)                 # only the transcript row is kept
  expect_equal(bed$V2, "100")                 # 1-based 101 -> 0-based 100
  expect_equal(bed$V3, "200")
  expect_equal(bed$V4, "ENST1")               # name taken from transcript_id
})

test_that(".normalize_transcripts errors on a GTF with no transcript rows", {
  gtf <- data.frame(seqnames = "chr1", start = 1, end = 9, strand = "+",
                    type = "exon", transcript_id = "ENST1")
  expect_error(.normalize_transcripts(gtf),
               class = "rnapeaks_error_not_found")
})

# --- intersect_peaks (overlap logic) --------------------------------------

test_that("fraction = 1 keeps only peaks fully contained on the same strand", {
  peaks <- rbind(bed6("chr1", 100, 200, "+"),   # inside the transcript
                 bed6("chr1", 500, 600, "-"))   # wrong strand
  tx    <- bed6("chr1", 0, 1000, "+")
  out <- suppressMessages(intersect_peaks(peaks, tx, fraction = 1,
                                          same_strand = TRUE))
  expect_equal(nrow(out), 1L)
  expect_equal(out$V2, "100")                 # the contained + strand peak
  expect_equal(out$V8, "0")                   # matched transcript start
})

test_that("a partial overlap is kept or dropped by the fraction threshold", {
  peaks <- bed6("chr1", 100, 200, "+")        # 1-based 101-200, width 100
  tx    <- bed6("chr1", 150, 1000, "+")       # overlap 151-200 = 50 bp -> 0.5
  expect_equal(nrow(suppressMessages(
    intersect_peaks(peaks, tx, fraction = 0.5))), 1L)
  expect_equal(nrow(suppressMessages(
    intersect_peaks(peaks, tx, fraction = 0.9))), 0L)
})

test_that("same_strand = FALSE lets an opposite-strand overlap through", {
  peaks <- bed6("chr1", 500, 600, "-")
  tx    <- bed6("chr1", 0, 1000, "+")
  expect_equal(nrow(suppressMessages(
    intersect_peaks(peaks, tx, fraction = 1, same_strand = FALSE))), 1L)
  expect_equal(nrow(suppressMessages(
    intersect_peaks(peaks, tx, fraction = 1, same_strand = TRUE))), 0L)
})

test_that("intersect_peaks aborts when peaks and transcripts share no chromosome", {
  expect_error(
    suppressMessages(intersect_peaks(bed6("chr1", 1, 9, "+"),
                                     bed6("chr2", 1, 9, "+"))),
    class = "rnapeaks_error_not_found"
  )
})

test_that("intersect_peaks validates required args and the fraction range", {
  expect_error(intersect_peaks(transcripts = bed6("chr1", 1, 9, "+")),
               class = "rnapeaks_error_invalid_arg")
  expect_error(intersect_peaks(bed6("chr1", 1, 9, "+")),
               class = "rnapeaks_error_invalid_arg")
  expect_error(intersect_peaks(bed6("chr1", 1, 9, "+"),
                               bed6("chr1", 1, 9, "+"), fraction = 2),
               class = "rnapeaks_error")
})
