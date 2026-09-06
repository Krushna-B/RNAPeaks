# Tests for plot_utr_binding() and its helpers in R/plot_utr_binding.R

# --- .prep_utr_bed_tracks -------------------------------------------------

test_that(".prep_utr_bed_tracks reduces a single data frame and names it", {
  bed <- make_raw_bed(n = 2, chr = "chr1", start = c(100, 120),
                      end = c(150, 200), strand = "+")
  res <- .prep_utr_bed_tracks(bed, default_name = "mybed")
  expect_named(res, "mybed")
  expect_s4_class(res$mybed, "GRanges")
  expect_equal(length(res$mybed), 1L)                       # overlapping peaks merged
  expect_equal(GenomicRanges::start(res$mybed), 100)
  expect_equal(GenomicRanges::end(res$mybed), 200)
})

test_that(".prep_utr_bed_tracks keeps supplied names for a list of tracks", {
  res <- .prep_utr_bed_tracks(list(A = make_raw_bed(1), B = make_raw_bed(1, chr = "chr2")))
  expect_named(res, c("A", "B"))
})

test_that(".prep_utr_bed_tracks rejects duplicate names and bad elements", {
  expect_error(.prep_utr_bed_tracks(list(A = make_raw_bed(1), A = make_raw_bed(1))),
               class = "rnapeaks_error_invalid_arg")
  expect_error(.prep_utr_bed_tracks(list(make_raw_bed(1), 42)),
               class = "rnapeaks_error_invalid_arg")
  expect_error(.prep_utr_bed_tracks(42), class = "rnapeaks_error_invalid_arg")
  expect_error(.prep_utr_bed_tracks(list()), class = "rnapeaks_error_invalid_arg")
})

# --- plot_utr_side_map ----------------------------------------------------

one_track_side <- function() {
  d <- data.frame(position_in_region = 1:100,
                  frequency = c(rep(0, 20), rep(1, 10), rep(0, 70)),
                  track = "A", gene_group = "All genes", n_events = 5L,
                  stringsAsFactors = FALSE)
  d$moving_avg <- d$frequency
  d
}

test_that("plot_utr_side_map returns a ggplot for a single track", {
  p <- plot_utr_side_map(one_track_side(), event_schema_utr, utr_style(),
                         "utr5", title = "5' UTR")
  expect_s3_class(p, "ggplot")
})

test_that("plot_utr_side_map aborts when a named palette misses a track", {
  d <- rbind(one_track_side(),
             transform(one_track_side(), track = "B"))
  style <- utr_style(palette = c(A = "red"))   # no entry for track B
  expect_error(plot_utr_side_map(d, event_schema_utr, style, "utr5"),
               class = "rnapeaks_error_invalid_arg")
})

test_that("plot_utr_side_map draws a group linetype legend for multiple groups", {
  d <- rbind(one_track_side(),
             transform(one_track_side(), gene_group = "Low"))
  d$gene_group[d$gene_group == "All genes"] <- "High"
  p <- plot_utr_side_map(d, event_schema_utr, utr_style(), "utr5",
                         show_group_legend = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that(".resolve_utr_group_spec rejects both args and bad group lists", {
  expect_error(.resolve_utr_group_spec(list(A = "X"), "Y"),
               class = "rnapeaks_error_invalid_arg")
  expect_error(.resolve_utr_group_spec(list("X"), NULL),          # unnamed
               class = "rnapeaks_error_invalid_arg")
  expect_error(.resolve_utr_group_spec(list(A = character(0)), NULL),
               class = "rnapeaks_error_invalid_arg")
  expect_null(.resolve_utr_group_spec(NULL, NULL))
  expect_equal(.resolve_utr_group_spec(NULL, "CXCR4"),
               list(`All genes` = "CXCR4"))
})

# --- gene_groups from a data frame / file ---------------------------------

test_that(".gene_group_table_to_list splits genes by group in first-seen order", {
  df <- data.frame(g = c("A", "B", "C", "D"),
                   grp = c("High", "High", "Low", "High"),
                   stringsAsFactors = FALSE)
  spec <- .gene_group_table_to_list(df)
  expect_equal(names(spec), c("High", "Low"))
  expect_equal(spec$High, c("A", "B", "D"))
  expect_equal(spec$Low, "C")
})

test_that(".gene_group_table_to_list drops blank rows and needs two columns", {
  df <- data.frame(g = c("A", "", "C"), grp = c("High", "Low", ""),
                   stringsAsFactors = FALSE)
  expect_equal(.gene_group_table_to_list(df), list(High = "A"))
  expect_error(.gene_group_table_to_list(data.frame(x = 1)),
               class = "rnapeaks_error_invalid_arg")
  expect_error(.gene_group_table_to_list(data.frame(g = "", grp = "")),
               class = "rnapeaks_error_invalid_arg")
})

test_that(".read_gene_groups_file auto-detects delimiter and header", {
  p1 <- tempfile(fileext = ".csv"); on.exit(unlink(p1), add = TRUE)
  writeLines(c("gene,group", "A,High", "B,Low"), p1)
  df1 <- .read_gene_groups_file(p1)
  expect_equal(nrow(df1), 2L)                       # header row dropped
  expect_equal(df1[[1]], c("A", "B"))
  expect_equal(df1[[2]], c("High", "Low"))

  p2 <- tempfile(fileext = ".tsv"); on.exit(unlink(p2), add = TRUE)
  writeLines(c("A\tHigh", "B\tLow"), p2)            # tab, no header
  df2 <- .read_gene_groups_file(p2)
  expect_equal(nrow(df2), 2L)
  expect_equal(df2[[1]], c("A", "B"))
})

test_that(".read_gene_groups_file errors on missing / empty / one-column files", {
  expect_error(.read_gene_groups_file("/no/such/file.csv"),
               class = "rnapeaks_error_not_found")
  p <- tempfile(); on.exit(unlink(p), add = TRUE)
  file.create(p)                                    # empty
  expect_error(.read_gene_groups_file(p), class = "rnapeaks_error_invalid_arg")
  writeLines(c("A", "B"), p)                        # single column
  expect_error(.read_gene_groups_file(p), class = "rnapeaks_error_invalid_arg")
})

test_that(".normalize_gene_groups routes list / data frame / path", {
  L <- list(High = "A", Low = "B")
  expect_identical(.normalize_gene_groups(L), L)    # list passes through
  df <- data.frame(g = c("A", "B"), grp = c("High", "Low"),
                   stringsAsFactors = FALSE)
  expect_equal(.normalize_gene_groups(df), list(High = "A", Low = "B"))
})

test_that(".resolve_utr_group_spec accepts a two-column data frame", {
  df <- data.frame(g = c("A", "B", "C"), grp = c("High", "Low", "High"),
                   stringsAsFactors = FALSE)
  expect_equal(.resolve_utr_group_spec(df, NULL),
               list(High = c("A", "C"), Low = "B"))
})

# --- plot_utr_binding entry point -----------------------------------------

test_that("plot_utr_binding reports missing bed / bad species via the error boundary", {
  err <- expect_error(suppressMessages(plot_utr_binding()))
  expect_match(conditionMessage(err), "Failed to generate UTR binding plot")
  expect_s3_class(err$parent, "rnapeaks_error_invalid_arg")

  err2 <- expect_error(
    suppressMessages(plot_utr_binding(make_raw_bed(), species = "nope"))
  )
  expect_s3_class(err2$parent, "rnapeaks_error_invalid_arg")
})

test_that("plot_utr_binding returns 5' and 3' plot/data pairs", {
  skip_if_not_installed("rtracklayer")
  gtf <- write_utr_gtf(); on.exit(unlink(gtf), add = TRUE)
  bed <- make_raw_bed(n = 2, chr = "chr1", start = c(110, 410),
                      end = c(130, 430), strand = "+", protein = "SRSF1")
  res <- suppressMessages(plot_utr_binding(bed, gtf = gtf))
  expect_s3_class(res$utr5$plot, "ggplot")
  expect_s3_class(res$utr3$plot, "ggplot")
  expect_s3_class(res$utr5$data, "data.frame")
  expect_equal(nrow(res$utr5$data), event_schema_utr$n_bins)
})
