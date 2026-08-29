# Tests for generate_control_peaks() and its helpers in
# R/generate_control_peaks.R
#
# The table readers and the peak-regrouping helper are deterministic and tested
# directly against hand-built inputs. The entry point is checked for argument
# validation, the no-shared-chromosome early return, and the overall 7-column
# contract of a full run (the C++ sampling engine may legitimately return zero
# controls for a tiny fixture, so the shape -- not the row count -- is asserted).

CONTROL_COLS <- c("chr", "peak_start", "peak_end", "name",
                  "strand", "control_start", "control_end")

# --- .empty_control_df ----------------------------------------------------

test_that(".empty_control_df has the final column shape and zero rows", {
  e <- .empty_control_df()
  expect_equal(nrow(e), 0L)
  expect_equal(names(e), CONTROL_COLS)
  expect_type(e$peak_start, "integer")
})

# --- read_annotation ------------------------------------------------------

test_that("read_annotation parses transcript/region from name when absent", {
  df <- data.frame(
    chr = "chr1", start = c(1000L, 2000L), end = c(1100L, 2200L),
    name = c("ENST1_CDS_1", "ENST1_UTR3_2"), strand = "+", gene = "G1",
    stringsAsFactors = FALSE
  )
  dt <- read_annotation(df)
  expect_equal(dt$transcript, c("ENST1", "ENST1"))
  expect_equal(dt$region, c("CDS", "UTR3"))
  expect_type(dt$start, "integer")
})

test_that("read_annotation reads a positional 7-column BED", {
  df <- stats::setNames(
    data.frame("chr1", 1000, 1100, "ENST2_CDS_1", ".", "+", "G1",
               stringsAsFactors = FALSE),
    paste0("c", 1:7)
  )
  dt <- read_annotation(df)
  expect_equal(dt$transcript, "ENST2")
  expect_equal(dt$region, "CDS")
  expect_equal(dt$gene, "G1")
})

test_that("read_annotation rejects too-few positional columns and bad names", {
  short <- stats::setNames(
    data.frame("chr1", 1, 9, "x", ".", "+", stringsAsFactors = FALSE),
    paste0("c", 1:6)
  )
  expect_error(read_annotation(short), class = "rnapeaks_error_invalid_bed")
  unparseable <- data.frame(chr = "chr1", start = 1L, end = 9L,
                            name = "NOSEP", strand = "+", gene = "G1")
  expect_error(read_annotation(unparseable),
               class = "rnapeaks_error_invalid_bed")
})

# --- read_genes -----------------------------------------------------------

test_that("read_genes reads named and positional gene tables", {
  named <- read_genes(data.frame(chr = "chr1", start = 0L, end = 5000L,
                                 gene = "G1"))
  expect_equal(as.character(named$gene), "G1")
  expect_type(named$start, "integer")

  positional <- read_genes(stats::setNames(
    data.frame("chr1", 0, 5000, "G1", stringsAsFactors = FALSE),
    paste0("c", 1:4)
  ))
  expect_equal(as.character(positional$gene), "G1")
})

test_that("read_genes rejects a positional table with fewer than four columns", {
  short <- stats::setNames(data.frame("chr1", 0, 5000, stringsAsFactors = FALSE),
                           paste0("c", 1:3))
  expect_error(read_genes(short), class = "rnapeaks_error_invalid_bed")
})

# --- .group_peaks_for_engine ----------------------------------------------

test_that(".group_peaks_for_engine collapses peaks and pools their transcripts", {
  # intersect_peaks output shape: V7 = fold change, V12 = transcript name.
  # Rows 1-2 are the same peak against two transcripts; row 3 is a second peak.
  df <- data.frame(
    V1 = "chr1", V2 = c(100, 100, 500), V3 = c(200, 200, 600),
    V4 = ".", V5 = ".", V6 = ".", V7 = c("2.5", "2.5", "1.1"),
    V8 = ".", V9 = ".", V10 = ".", V11 = ".",
    V12 = c("ENST1", "ENST2", "ENST3"), V13 = ".", V14 = ".",
    stringsAsFactors = FALSE
  )
  grouped <- .group_peaks_for_engine(df)
  expect_equal(nrow(grouped), 2L)
  expect_equal(grouped$peak_range, c("100-200", "500-600"))   # original order
  expect_equal(grouped$start, c(100L, 500L))
  expect_equal(grouped$FC, c("2.5", "1.1"))
  expect_equal(grouped$transcripts[[1]], c("ENST1", "ENST2")) # pooled per peak
  expect_equal(grouped$transcripts[[2]], "ENST3")
})

# --- generate_control_peaks -----------------------------------------------

# BED8 peak (chr, start, end, name, gene, strand, FC, pval) fully inside ENST1.
peaks8 <- function() data.frame(
  V1 = "chr1", V2 = 1000, V3 = 1100, V4 = "pk", V5 = "G1",
  V6 = "+", V7 = "5.0", V8 = "0.01", stringsAsFactors = FALSE
)
tx6 <- function() data.frame(
  V1 = "chr1", V2 = 0, V3 = 5000, V4 = "ENST1", V5 = ".", V6 = "+",
  stringsAsFactors = FALSE
)
anno_df <- function(chr = "chr1") data.frame(
  chr = chr, start = c(1000L, 2000L, 3000L, 4000L),
  end = c(1100L, 2500L, 3500L, 4800L),
  name = paste0("ENST1_CDS_", 1:4), strand = "+", gene = "G1",
  stringsAsFactors = FALSE
)
gene_df <- function() data.frame(chr = "chr1", start = 0L, end = 5000L,
                                 gene = "G1", stringsAsFactors = FALSE)

test_that("generate_control_peaks validates required arguments and threads", {
  expect_error(generate_control_peaks(anno = anno_df(), gene = gene_df(),
                                      transcripts = tx6()),
               class = "rnapeaks_error_invalid_arg")               # raw_peaks
  expect_error(generate_control_peaks(peaks8(), anno_df(), gene_df(), tx6(),
                                      threads = 0),
               class = "rnapeaks_error_invalid_arg")
  expect_error(generate_control_peaks(42, anno_df(), gene_df(), tx6()),
               class = "rnapeaks_error_invalid_arg")               # wrong type
})

test_that("generate_control_peaks returns an empty result off a shared chromosome", {
  # Peaks/transcripts on chr1, annotation on chr2 -> no common chromosome.
  res <- suppressMessages(
    generate_control_peaks(peaks8(), anno_df("chr2"), gene_df(), tx6())
  )
  expect_equal(nrow(res), 0L)
  expect_equal(names(res), CONTROL_COLS)
})

test_that("generate_control_peaks runs the full pipeline to the 7-column contract", {
  res <- suppressMessages(
    generate_control_peaks(peaks8(), anno_df(), gene_df(), tx6(), seed = 42)
  )
  expect_equal(names(res), CONTROL_COLS)
})
