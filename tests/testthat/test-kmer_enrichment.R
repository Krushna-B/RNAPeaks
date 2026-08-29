# Tests for kmer_enrichment() and its helpers in R/kmer_enrichment.R
#
# The counting / table / plotting helpers are genome-free and hand-derivable and
# are tested directly. .ids_to_span_granges resolves ids against the shared
# make_gtf() fixture (spans are hand-derived from those coordinates). The full
# entry point and .resolve_kmer_set need a BSgenome, so they are gated.

skip_no_hg38 <- function() skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

kbed <- function(chr, start, end, strand = "+") {
  data.frame(chr = chr, start = start, end = end, name = "p", score = 0L,
             strand = strand, stringsAsFactors = FALSE)
}

# --- .is_bed_set ----------------------------------------------------------

test_that(".is_bed_set recognizes BED-shaped inputs, not id vectors", {
  expect_true(.is_bed_set(kbed("chr1", 1, 9)))
  expect_true(.is_bed_set(list(kbed("chr1", 1, 9), kbed("chr2", 1, 9))))
  expect_false(.is_bed_set(c("ENSG1", "ENSG2")))       # ids, not a path
  expect_false(.is_bed_set("/no/such/file.bed"))
})

# --- .kmer_freq -----------------------------------------------------------

test_that(".kmer_freq counts overlapping k-mers and normalizes to 1", {
  # "AAAA" contains AA at positions 1-2, 2-3, 3-4 -> 3 counts, all mass on AA.
  kf <- .kmer_freq(Biostrings::DNAStringSet("AAAA"), k = 2)
  expect_equal(unname(kf$count["AA"]), 3)
  expect_equal(unname(kf$freq["AA"]), 1)
  expect_equal(sum(kf$freq), 1)
})

test_that(".kmer_freq errors on empty input or when no valid k-mer exists", {
  expect_error(.kmer_freq(Biostrings::DNAStringSet(), k = 2),
               class = "rnapeaks_error_not_found")
  # a single base is shorter than k -> no 2-mers -> zero total
  expect_error(.kmer_freq(Biostrings::DNAStringSet("N"), k = 2),
               class = "rnapeaks_error_not_found")
})

# --- .build_kmer_table ----------------------------------------------------

test_that(".build_kmer_table joins, differences, and ranks by enrichment", {
  kf_a <- list(count = c(AA = 3, AC = 1), freq = c(AA = 0.75, AC = 0.25))
  kf_b <- list(count = c(AA = 1, AC = 3), freq = c(AA = 0.25, AC = 0.75))
  tbl <- .build_kmer_table(kf_a, kf_b)
  # AA is enriched in A (+0.5), AC in B (-0.5) -> AA ranks first.
  expect_equal(tbl$kmer, c("AA", "AC"))
  expect_equal(tbl$difference, c(0.5, -0.5))
  expect_equal(tbl$freq_b, c(0.25, 0.75))
  expect_equal(tbl$rank, c(1L, 2L))
})

# --- plots ----------------------------------------------------------------

test_that("the k-mer plotters return ggplots for labelled and unlabelled cases", {
  tbl <- .build_kmer_table(
    list(count = c(AA = 3, AC = 1), freq = c(AA = 0.75, AC = 0.25)),
    list(count = c(AA = 1, AC = 3), freq = c(AA = 0.25, AC = 0.75))
  )
  expect_s3_class(.plot_kmer_scatter(tbl, "A", "B", top_n = 1, ""), "ggplot")
  expect_s3_class(.plot_kmer_scatter(tbl, "A", "B", top_n = 0, ""), "ggplot")
  expect_s3_class(.plot_kmer_rank(tbl, "A", "B", ""), "ggplot")
})

# --- .bed_set_to_granges --------------------------------------------------

test_that(".bed_set_to_granges reduces a single BED and a list of BEDs", {
  one <- .bed_set_to_granges(kbed("chr1", 100, 200))
  expect_s4_class(one, "GRanges")
  expect_length(one, 1L)
  many <- .bed_set_to_granges(list(kbed("chr1", 100, 200),
                                   kbed("chr1", 500, 600)))
  expect_length(many, 2L)                    # two disjoint ranges, not merged
})

# --- .ids_to_span_granges -------------------------------------------------

test_that(".ids_to_span_granges spans each id from its min start to max end", {
  gtf <- make_gtf()
  # ENST00000003 (U2AF2): 5000-5600 (+); ENSG00000136450 (SRSF1) longest tx
  # ENST00000001: 100-1100 (+).
  gr <- .ids_to_span_granges(c("ENST00000003", "ENSG00000136450"), gtf)
  expect_equal(GenomicRanges::start(gr), c(5000, 100))
  expect_equal(GenomicRanges::end(gr),   c(5600, 1100))
  expect_equal(as.character(GenomicRanges::strand(gr)), c("+", "+"))
})

# --- kmer_enrichment (validation before the genome is resolved) -----------

test_that("kmer_enrichment reports missing set / k arguments", {
  expect_error(kmer_enrichment(set_b = kbed("chr1", 1, 9), k = 2),
               class = "rnapeaks_error_invalid_arg")
  expect_error(kmer_enrichment(kbed("chr1", 1, 9), kbed("chr1", 1, 9)),
               class = "rnapeaks_error_invalid_arg")   # k missing
})

test_that("kmer_enrichment caps k at 12", {
  expect_error(
    kmer_enrichment(kbed("chr1", 1, 9), kbed("chr1", 1, 9), k = 13),
    class = "rnapeaks_error_invalid_arg"
  )
})

# --- kmer_enrichment (full run, needs a genome) ---------------------------

test_that("kmer_enrichment returns a ranked table and two plots for BED sets", {
  skip_no_hg38()
  res <- suppressMessages(kmer_enrichment(
    kbed("chr1", 100000, 100200), kbed("chr1", 200000, 200200), k = 2
  ))
  expect_named(res, c("table", "plots"))
  expect_true(all(c("kmer", "freq_a", "freq_b", "difference", "rank") %in%
                    names(res$table)))
  expect_equal(res$table$rank, seq_len(nrow(res$table)))   # ranked 1..n
  expect_s3_class(res$plots$scatter, "ggplot")
  expect_s3_class(res$plots$rank, "ggplot")
})
