# Tests for the event schemas in R/sm-event_schema.R
#
# build_regions() places a window of width_exon into each exon and width_intron
# into each intron around the event's splice junctions. Coordinates below are
# hand-derived for width_exon = 10, width_intron = 20 (region_width = 31).

WE <- 10L; WI <- 20L

# --- Skipped Exon (4 regions) ---------------------------------------------
# upstream exon 1000-1100, skipped exon 1400-1500, downstream exon 1800-1900
se_event <- function(strand = "+") {
  data.frame(chr = "1", strand = strand,
             upstreamES = 1000, upstreamEE = 1100,
             exonStart_0base = 1400, exonEnd = 1500,
             downstreamES = 1800, downstreamEE = 1900,
             stringsAsFactors = FALSE)
}

test_that("SE build_regions makes 4 windows straddling the four junctions", {
  gr <- event_schema_se$build_regions(se_event(), WE, WI)
  expect_length(gr, 4)
  expect_equal(GenomicRanges::start(gr), c(1090, 1380, 1490, 1780))
  expect_equal(GenomicRanges::end(gr),   c(1120, 1410, 1520, 1810))
  expect_equal(GenomicRanges::mcols(gr)$region_idx, 1:4)
  expect_true(all(GenomicRanges::mcols(gr)$event_id == 1L))
})

test_that("SE region widths equal width_exon + width_intron + 1", {
  expect_equal(event_schema_se$region_width(WE, WI), 31L)
  gr <- event_schema_se$build_regions(se_event(), WE, WI)
  expect_true(all(IRanges::width(gr) == 31L))
})

test_that("SE build_regions is strand-agnostic (reorientation happens in the scorer)", {
  plus  <- event_schema_se$build_regions(se_event("+"), WE, WI)
  minus <- event_schema_se$build_regions(se_event("-"), WE, WI)
  expect_equal(GenomicRanges::start(plus), GenomicRanges::start(minus))
  expect_equal(as.character(GenomicRanges::strand(minus)), rep("-", 4))
})

test_that("SE build_regions returns empty GRanges for no events", {
  expect_length(event_schema_se$build_regions(se_event()[0, ], WE, WI), 0)
})

test_that("SE plot_layout spaces regions with the fixed gap", {
  lay <- event_schema_se$plot_layout(WE, WI)     # bin_width 30, gap 80
  expect_equal(lay$region_starts, c(0, 110, 140, 250))
  expect_equal(lay$x_max, 280)
  expect_equal(lay$boundary_offsets, c(10, 20, 10, 20))
})

# --- Retained Intron (2 regions) ------------------------------------------
# upstream exon 1000-1100, downstream exon 1500-1600 (intron between)
test_that("RI build_regions makes 2 windows around the retained intron", {
  ev <- data.frame(chr = "1", strand = "+",
                   upstreamES = 1000, upstreamEE = 1100,
                   downstreamES = 1500, downstreamEE = 1600)
  gr <- event_schema_ri$build_regions(ev, WE, WI)
  expect_length(gr, 2)
  expect_equal(GenomicRanges::start(gr), c(1090, 1480))
  expect_equal(GenomicRanges::end(gr),   c(1120, 1510))
})

# --- Alternative 5' splice site (+ strand) --------------------------------
# alt (short) exon 1000-1100, flanking exon 1500-1600
test_that("a5ss build_regions anchors on the short donor and flanking exon (+)", {
  ev <- data.frame(chr = "1", strand = "+",
                   longExonStart_0base = 1000, longExonEnd = 1200,
                   shortES = 1000, shortEE = 1100,
                   flankingES = 1500, flankingEE = 1600)
  gr <- event_schema_a5ss$build_regions(ev, WE, WI)
  expect_length(gr, 2)
  expect_equal(GenomicRanges::start(gr), c(1090, 1480))
  expect_equal(GenomicRanges::end(gr),   c(1120, 1510))
})

# --- Alternative 3' splice site (+ strand) --------------------------------
# flanking exon 1000-1100, alt (short) exon 1500-1600
test_that("a3ss build_regions anchors on the flanking exon and short acceptor (+)", {
  ev <- data.frame(chr = "1", strand = "+",
                   longExonStart_0base = 1400, longExonEnd = 1600,
                   shortES = 1500, shortEE = 1600,
                   flankingES = 1000, flankingEE = 1100)
  gr <- event_schema_a3ss$build_regions(ev, WE, WI)
  expect_length(gr, 2)
  expect_equal(GenomicRanges::start(gr), c(1090, 1480))
  expect_equal(GenomicRanges::end(gr),   c(1120, 1510))
})

# --- schema metadata ------------------------------------------------------

test_that("all trio schemas expose the expected structure", {
  for (sch in list(event_schema_se, event_schema_ri,
                   event_schema_a5ss, event_schema_a3ss)) {
    expect_true(is.function(sch$build_regions))
    expect_equal(sch$group_set, c("Negative", "Positive", "Control"))
    expect_true(all(c("chr", "strand", "GeneID", "PValue", "FDR",
                      "IncLevelDifference") %in% sch$required_cols))
  }
  expect_equal(event_schema_se$n_regions, 4L)
  expect_equal(event_schema_ri$n_regions, 2L)
})

# --- .schematic_ext -------------------------------------------------------
# The a5ss / a3ss cartoon draws the long-isoform extension. Its length is the
# rounded median real extension, capped to the intron side of the window; with
# no usable stat it falls back to max(8, 12% of bin_width).
# a5ss layout at (10, 20): boundary_offsets = c(10, 20) so intron_w = 20,
# bin_width = 30 -> fallback = max(8, round(3.6)) = 8.

test_that(".schematic_ext uses the rounded median extension when available", {
  layout <- event_schema_a5ss$plot_layout(10, 20)
  expect_equal(.schematic_ext(list(median_ext = 15.4), layout), 15)   # rounds
})

test_that(".schematic_ext caps the extension at the intron width", {
  layout <- event_schema_a5ss$plot_layout(10, 20)
  expect_equal(.schematic_ext(list(median_ext = 50), layout), 20)     # capped
})

test_that(".schematic_ext falls back when the stat is missing or non-positive", {
  small <- event_schema_a5ss$plot_layout(10, 20)     # bin_width 30 -> floor 8
  large <- event_schema_a5ss$plot_layout(50, 250)    # bin_width 300 -> 12% = 36
  expect_equal(.schematic_ext(NULL, small), 8L)
  expect_equal(.schematic_ext(list(median_ext = 0), small), 8L)
  expect_equal(.schematic_ext(list(median_ext = NA_real_), large), 36)
})
