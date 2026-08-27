# Tests for build_utr_events() + resolve_transcript_ids() in R/utr-build_events.R
#
# Verifies the strand-aware 5'/3' assignment, the longest-per-(gene, side)
# transcript pick, and the piece tables -- using make_utr_gtf() whose expected
# outcomes are documented in the fixture.

# --- side assignment (+ strand) -------------------------------------------

test_that("on + strand, UTR upstream of the CDS is 5' and downstream is 3'", {
  built <- build_utr_events(make_utr_gtf())
  g1 <- built$events[built$events$gene_id == "ENSG1", ]
  expect_equal(g1$tx5, "ENST1")
  expect_equal(g1$tx3, "ENST1")
  expect_equal(g1$utr5_len, 51)   # 100-150
  expect_equal(g1$utr3_len, 51)   # 400-450

  p5 <- built$utr5_pieces[built$utr5_pieces$event_idx == which(built$events$gene_id == "ENSG1"), ]
  expect_equal(p5$start, 100)
  expect_equal(p5$end,   150)
})

# --- side assignment (- strand) -------------------------------------------

test_that("on - strand the 5'/3' sides swap relative to genomic coordinates", {
  built <- build_utr_events(make_utr_gtf())
  idx <- which(built$events$gene_id == "ENSG2")
  # downstream-of-CDS piece (400-450) is the 5' UTR on the minus strand
  p5 <- built$utr5_pieces[built$utr5_pieces$event_idx == idx, ]
  p3 <- built$utr3_pieces[built$utr3_pieces$event_idx == idx, ]
  expect_equal(p5$start, 400)
  expect_equal(p3$start, 100)
})

# --- longest transcript per (gene, side) ----------------------------------

test_that("each side independently picks the transcript with the longest UTR", {
  built <- build_utr_events(make_utr_gtf())
  g3 <- built$events[built$events$gene_id == "ENSG3", ]
  expect_equal(g3$tx5, "ENST3")   # 5' UTR 51 bp beats ENST4's 21 bp
  expect_equal(g3$tx3, "ENST4")   # 3' UTR 61 bp beats ENST3's 51 bp
  expect_equal(g3$utr5_len, 51)
  expect_equal(g3$utr3_len, 61)
})

test_that("only the chosen transcript's pieces appear in the piece tables", {
  built <- build_utr_events(make_utr_gtf())
  idx <- which(built$events$gene_id == "ENSG3")
  p5 <- built$utr5_pieces[built$utr5_pieces$event_idx == idx, ]
  # ENST3 was chosen for 5'; its UTR starts at 100 (ENST4's would start at 130)
  expect_equal(p5$start, 100)
})

# --- error paths ----------------------------------------------------------

test_that("build_utr_events validates the transcripts filter", {
  expect_error(build_utr_events(make_utr_gtf(), transcripts = 1L),
               class = "rnapeaks_error_invalid_arg")
})

test_that("build_utr_events aborts when no CDS/UTR rows survive", {
  gtf <- make_utr_gtf(); gtf$type <- "exon"
  expect_error(build_utr_events(gtf), class = "rnapeaks_error_not_found")
})

test_that("build_utr_events aborts when no UTR shares a transcript with a CDS", {
  gtf <- make_utr_gtf()[make_utr_gtf()$type == "UTR", ]  # UTRs but no CDS
  expect_error(build_utr_events(gtf), class = "rnapeaks_error_not_found")
})

# --- resolve_transcript_ids -----------------------------------------------

test_that("resolve_transcript_ids accepts transcript ids, gene ids, and symbols", {
  gtf <- make_utr_gtf()
  expect_equal(resolve_transcript_ids(gtf, "ENST1"), "ENST1")
  expect_setequal(resolve_transcript_ids(gtf, "ENSG3"), c("ENST3", "ENST4"))
  expect_equal(resolve_transcript_ids(gtf, "GENEA"), "ENST1")
})

test_that("resolve_transcript_ids tolerates some unmatched ids but aborts on none", {
  gtf <- make_utr_gtf()
  expect_equal(suppressMessages(resolve_transcript_ids(gtf, c("ENST1", "NOPE"))), "ENST1")
  expect_error(resolve_transcript_ids(gtf, "NOPE"), class = "rnapeaks_error_not_found")
})

test_that("build_utr_events honours the transcripts filter", {
  built <- build_utr_events(make_utr_gtf(), transcripts = "ENST1")
  expect_equal(built$events$gene_id, "ENSG1")
})
