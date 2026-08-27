# Tests for select_transcript() + its helpers detect_species() / normalize_label()
# in R/gene-select.R
#
# Contract:
#   - detect_species: "mouse" if any gene_id starts ENSMUSG, else "human" if any
#     starts ENSG, else abort_species_unknown.
#   - normalize_label: human -> all caps; mouse -> Title case.
#   - select_transcript: filter by gene (symbol or Ensembl id), then pick a
#     transcript (explicit TxID, or the LONGEST transcript by width).

# --- detect_species -------------------------------------------------------

test_that("detect_species reads the species from gene_id prefixes", {
  expect_equal(detect_species(make_gtf()), "human")
  expect_equal(detect_species(data.frame(gene_id = "ENSMUSG000001")), "mouse")
})

test_that("detect_species aborts when no gene_id prefix is recognised", {
  expect_error(detect_species(data.frame(gene_id = "FBgn0000001")),
               class = "rnapeaks_error_species_unknown")
})

# --- normalize_label ------------------------------------------------------

test_that("normalize_label upper-cases for human and Title-cases for mouse", {
  expect_equal(normalize_label("srsf1", "human"), "SRSF1")
  expect_equal(normalize_label("SRSF1", "mouse"), "Srsf1")
  expect_equal(normalize_label("myGene", "mouse"), "Mygene")
})

# --- select_transcript: validation ----------------------------------------

test_that("select_transcript validates gtf shape and required columns", {
  expect_error(select_transcript(data.frame()), class = "rnapeaks_error_invalid_arg")
  bad <- make_gtf()[, setdiff(names(make_gtf()), "width")]
  expect_error(select_transcript(bad, geneID = "SRSF1"),
               class = "rnapeaks_error_invalid_gtf")
})

test_that("select_transcript validates geneID / TxID and requires at least one", {
  gtf <- make_gtf()
  expect_error(select_transcript(gtf, geneID = c("A", "B")), class = "rnapeaks_error_invalid_arg")
  expect_error(select_transcript(gtf, geneID = ""), class = "rnapeaks_error_invalid_arg")
  expect_error(select_transcript(gtf, geneID = 123), class = "rnapeaks_error_invalid_arg")
  expect_error(select_transcript(gtf, TxID = c("A", "B")), class = "rnapeaks_error_invalid_arg")
  expect_error(select_transcript(gtf), class = "rnapeaks_error_invalid_arg")  # neither given
})

# --- select_transcript: gene lookup + longest-transcript pick -------------

test_that("selecting by gene symbol returns the longest transcript", {
  out <- select_transcript(make_gtf(), geneID = "SRSF1")
  expect_equal(unique(out$transcript_id), "ENST00000001")  # width 1000 > 400
})

test_that("gene symbol match is case-insensitive and Ensembl gene ids also work", {
  expect_equal(unique(select_transcript(make_gtf(), geneID = "srsf1")$transcript_id),
               "ENST00000001")
  expect_equal(unique(select_transcript(make_gtf(), geneID = "ENSG00000136450")$transcript_id),
               "ENST00000001")
})

test_that("an unknown gene aborts as not-found", {
  expect_error(select_transcript(make_gtf(), geneID = "NOPE"),
               class = "rnapeaks_error_not_found")
})

# --- select_transcript: explicit transcript selection ---------------------

test_that("an explicit Ensembl TxID overrides the longest-transcript default", {
  out <- select_transcript(make_gtf(), TxID = "ENST00000002")  # the shorter one
  expect_equal(unique(out$transcript_id), "ENST00000002")
})

test_that("TxID lookup is case-insensitive for Ensembl ids", {
  out <- select_transcript(make_gtf(), TxID = "enst00000002")
  expect_equal(unique(out$transcript_id), "ENST00000002")
})

test_that("a transcript can be chosen by transcript_name", {
  out <- select_transcript(make_gtf(), TxID = "SRSF1-202")
  expect_equal(unique(out$transcript_id), "ENST00000002")
})

test_that("matching by transcript name without a transcript_name column aborts", {
  gtf <- make_gtf()[, setdiff(names(make_gtf()), "transcript_name")]
  expect_error(select_transcript(gtf, TxID = "SRSF1-202"),
               class = "rnapeaks_error_invalid_gtf")
})

test_that("an unknown transcript aborts as not-found", {
  expect_error(select_transcript(make_gtf(), TxID = "ENST09999999"),
               class = "rnapeaks_error_not_found")
})

test_that("geneID and TxID together narrow to the named transcript within the gene", {
  out <- select_transcript(make_gtf(), geneID = "SRSF1", TxID = "ENST00000002")
  expect_equal(unique(out$transcript_id), "ENST00000002")
})
