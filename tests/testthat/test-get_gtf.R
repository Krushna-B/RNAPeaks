# Tests for get_GTF() / normalize_gtf() / verify_gtf() in R/gene-get_gtf.R
#
# Contract:
#   - normalize_gtf: chr-normalizes seqnames, coerces type to character, and
#     collapses Ensembl five_prime_utr / three_prime_utr into a single "UTR".
#   - verify_gtf: passes (invisibly TRUE) for a non-empty frame with all
#     required columns; otherwise aborts with class rnapeaks_error_invalid_gtf.
#   - get_GTF: with `file` it validates/imports a path; otherwise it loads the
#     bundled dataset for a supported species.

# --- normalize_gtf --------------------------------------------------------

test_that("normalize_gtf canonicalises seqnames", {
  gtf <- data.frame(seqnames = c("chr1", "chrX", " chr2 "),
                    type = "exon", stringsAsFactors = FALSE)
  out <- normalize_gtf(gtf)
  expect_equal(out$seqnames, c("1", "X", "2"))
})

test_that("normalize_gtf collapses UTR feature types and coerces type to character", {
  gtf <- data.frame(
    seqnames = "1",
    type = factor(c("five_prime_utr", "three_prime_utr", "UTR", "exon")),
    stringsAsFactors = FALSE
  )
  out <- normalize_gtf(gtf)
  expect_type(out$type, "character")
  expect_equal(out$type, c("UTR", "UTR", "UTR", "exon"))
})

# --- verify_gtf -----------------------------------------------------------

test_that("verify_gtf passes for a well-formed GTF and returns invisibly", {
  expect_true(verify_gtf(make_gtf()))
  expect_invisible(verify_gtf(make_gtf()))
})

test_that("verify_gtf rejects non-frames, missing columns, and empty frames", {
  expect_error(verify_gtf(list(a = 1)), class = "rnapeaks_error_invalid_gtf")
  expect_error(verify_gtf(make_gtf()[, setdiff(names(make_gtf()), "gene_name")]),
               class = "rnapeaks_error_invalid_gtf")
  expect_error(verify_gtf(make_gtf()[0, ]), class = "rnapeaks_error_invalid_gtf")
})

# --- get_GTF: file argument validation ------------------------------------

test_that("get_GTF validates the file argument before touching disk", {
  expect_error(get_GTF("hg38", file = 123), class = "rnapeaks_error_invalid_arg")
  expect_error(get_GTF("hg38", file = NA), class = "rnapeaks_error_invalid_arg")
  expect_error(get_GTF("hg38", file = ""), class = "rnapeaks_error_invalid_arg")
})

test_that("get_GTF reports a missing file path as not-found", {
  expect_error(get_GTF("hg38", file = tempfile(fileext = ".gtf")),
               class = "rnapeaks_error_not_found")
})

# --- get_GTF: species argument validation ---------------------------------

test_that("get_GTF validates the species argument", {
  expect_error(get_GTF(123, file = NULL), class = "rnapeaks_error_invalid_arg")
  expect_error(get_GTF("hg99", file = NULL), class = "rnapeaks_error_invalid_arg")
})

# --- get_GTF: happy paths (heavier / dependency-guarded) ------------------

test_that("get_GTF imports a local GTF file and returns a verified, normalized frame", {
  skip_if_not_installed("rtracklayer")
  path <- tempfile(fileext = ".gtf")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    '1\tsrc\ttranscript\t100\t1100\t.\t+\t.\tgene_id "ENSG1"; gene_name "AAA"; transcript_id "ENST1";',
    '1\tsrc\texon\t100\t300\t.\t+\t.\tgene_id "ENSG1"; gene_name "AAA"; transcript_id "ENST1";'
  ), path)
  gtf <- get_GTF("hg38", file = path)
  expect_s3_class(gtf, "data.frame")
  expect_true(all(c("seqnames", "start", "end", "strand", "type",
                    "gene_id", "gene_name", "transcript_id") %in% names(gtf)))
  expect_equal(unique(as.character(gtf$seqnames)), "1")  # already normalized
})

test_that("get_GTF loads the bundled species dataset (normalized seqnames)", {
  skip_on_cran()
  skip_if_not(exists("gtf_hg38") || requireNamespace("RNAPeaks", quietly = TRUE),
              "bundled gtf_hg38 unavailable")
  gtf <- get_GTF("hg38", file = NULL)
  expect_s3_class(gtf, "data.frame")
  expect_gt(nrow(gtf), 0)
  # canonical seqnames carry no "chr" prefix
  expect_false(any(grepl("^chr", as.character(gtf$seqnames), ignore.case = TRUE)))
})
