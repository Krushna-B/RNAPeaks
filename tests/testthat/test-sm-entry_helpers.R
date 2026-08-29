# Tests for the shared entry-point helpers in R/sm-entry_helpers.R
#
# Covered here: the deterministic helpers that need no genome --
#   .peaks_to_granges, .resolve_genome, validate_sm_inputs, wrap_sm_errors,
#   .normalize_motifs, and the empty-input guard of .extract_region_seqs.
# The genome-backed helpers (.extract_region_seqs full path,
# .prepare_sequence_map_prep, .run_sequence_map) are exercised end to end in the
# sequence-map entry-point tests, gated on the BSgenome packages.

# --- .peaks_to_granges ----------------------------------------------------

test_that(".peaks_to_granges validates, 0-base-shifts, and reduces a BED", {
  # Two overlapping + strand peaks. check_bed treats starts as 0-based, so
  # start 100 -> 101; the overlap collapses to a single reduced range 101-250.
  raw <- data.frame(
    chr = "chr1", start = c(100, 150), end = c(200, 250),
    name = "P", score = 0L, strand = "+",
    stringsAsFactors = FALSE
  )
  gr <- .peaks_to_granges(raw)
  expect_s4_class(gr, "GRanges")
  expect_length(gr, 1L)
  expect_equal(GenomicRanges::start(gr), 101)
  expect_equal(GenomicRanges::end(gr), 250)
  expect_equal(as.character(GenomeInfoDb::seqnames(gr)), "1")   # chr normalized
})

test_that(".peaks_to_granges aborts when a BED path does not exist", {
  expect_error(.peaks_to_granges("/no/such/file.bed"),
               class = "rnapeaks_error_not_found")
})

# --- .resolve_genome ------------------------------------------------------

test_that(".resolve_genome rejects an unknown genome key", {
  expect_error(.resolve_genome("hg99"), class = "rnapeaks_error_invalid_arg")
})

test_that(".resolve_genome resolves a known key or reports the missing package", {
  # Whichever way the test machine is set up, exactly one branch is correct:
  # the BSgenome is loaded, or its absence is reported as a not-found error.
  pkg <- "BSgenome.Hsapiens.UCSC.hg38"
  if (requireNamespace(pkg, quietly = TRUE)) {
    expect_s4_class(.resolve_genome("hg38"), "BSgenome")
  } else {
    expect_error(.resolve_genome("hg38"), class = "rnapeaks_error_not_found")
  }
})

# --- validate_sm_inputs ---------------------------------------------------

test_that("validate_sm_inputs passes a fully valid argument set", {
  expect_silent(
    validate_sm_inputs(
      events = data.frame(), opts = splicing_options(),
      style = splicing_style(), title = "t",
      bed_file = data.frame(), sequence = "ACGT",
      genome = "hg38", motif_mode = "combined"
    )
  )
})

test_that("validate_sm_inputs rejects each malformed argument", {
  good <- list(events = data.frame(), opts = splicing_options(),
               style = splicing_style(), title = "t")
  # Direct assignment (not modifyList, which would recurse into list-valued args).
  err <- function(field, value) {
    args <- good
    args[[field]] <- value
    expect_error(do.call(validate_sm_inputs, args),
                 class = "rnapeaks_error_invalid_arg")
  }
  err("title",      5)                         # not a string
  err("events",     list())                    # not a data frame
  err("opts",       list())                    # not a splicing_options result
  err("style",      list())                    # not a splicing_style result
  err("bed_file",   42)                        # not path / df / list of df
  err("sequence",   c("ACG", ""))              # empty motif string
  err("genome",     42)                        # not NULL / string / BSgenome
  err("motif_mode", "sideways")                # not combined / individual
})

# --- wrap_sm_errors -------------------------------------------------------

test_that("wrap_sm_errors returns the body value on success", {
  expect_equal(wrap_sm_errors("splicing map", 40 + 2), 42)
})

test_that("wrap_sm_errors re-raises rnapeaks errors with just the context line", {
  # A classed rnapeaks error keeps its own message and gains only the context;
  # it must NOT be labelled an unexpected error.
  cnd <- tryCatch(wrap_sm_errors("splicing map", abort_invalid_arg("bad peaks")),
                  error = function(e) e)
  msg <- conditionMessage(cnd)
  expect_match(msg, "Failed to generate splicing map")
  expect_false(grepl("unexpected error", msg))
})

test_that("wrap_sm_errors annotates unexpected (non-rnapeaks) errors", {
  expect_error(wrap_sm_errors("sequence map", stop("boom")),
               regexp = "unexpected error")
})

# --- .extract_region_seqs (empty guard) -----------------------------------

test_that(".extract_region_seqs short-circuits on empty regions", {
  res <- .extract_region_seqs(GenomicRanges::GRanges(), genome = NULL)
  expect_length(res$regions, 0L)
  expect_s4_class(res$seqs, "DNAStringSet")
  expect_length(res$seqs, 0L)
})

# --- .normalize_motifs ----------------------------------------------------

test_that(".normalize_motifs upper-cases, trims, and converts U to T", {
  expect_message(out <- .normalize_motifs("acgu"), regexp = "Converting U")
  expect_equal(out, "ACGT")
  expect_equal(suppressMessages(.normalize_motifs(c(" AuG ", "cgn"))),
               c("ATG", "CGN"))
})

test_that(".normalize_motifs accepts IUPAC ambiguity codes without a U message", {
  expect_silent(out <- .normalize_motifs("RYSWKMBDHVN"))
  expect_equal(out, "RYSWKMBDHVN")
})

test_that(".normalize_motifs rejects non-IUPAC characters", {
  expect_error(.normalize_motifs("ACGX"), class = "rnapeaks_error_invalid_arg")
})
