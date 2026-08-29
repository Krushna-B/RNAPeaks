# Tests for the shared event-map pipeline in R/sm-event_map_pipeline.R
#
# The pure helpers (.assemble_frequency_frame, .smooth_by_group, ...) are
# checked with hand-built per_group inputs and hand-derived expected frames.
# prepare_event_map() and event_map_pipeline() are exercised end to end with a
# small SE fixture and a stub scorer / plot_fn, so the pipeline wiring is tested
# without depending on real genome scoring (covered in test-sm-scorers.R).

# A per_group list shaped exactly like the pipeline builds internally:
# named groups, each list(counts, n[, hits]).
per_group_fixture <- function() {
  list(
    Positive = list(counts = c(2, 0, 0, 4, 0, 0), n = 4L, hits = NULL),
    Negative = list(counts = c(3, 0, 0, 0, 0, 0), n = 5L, hits = NULL),
    Control  = list(counts = c(1, 1, 1, 1, 1, 1), n = 3L, hits = NULL)
  )
}

# --- .check_event_column_types --------------------------------------------

test_that("numeric event columns must be numeric; text columns are exempt", {
  ok <- data.frame(
    chr = "1", strand = "+", GeneID = "G1",
    PValue = 0.01, FDR = 0.01, IncLevelDifference = 0.3,
    exonStart_0base = 1400, exonEnd = 1500,
    upstreamES = 1000, upstreamEE = 1100,
    downstreamES = 1800, downstreamEE = 1900,
    stringsAsFactors = FALSE
  )
  expect_silent(.check_event_column_types(ok, event_schema_se))

  bad <- ok
  bad$PValue <- "0.01"   # a character coordinate/stat is the classic reader bug
  expect_error(.check_event_column_types(bad, event_schema_se),
               class = "rnapeaks_error_invalid_arg")
})

# --- .validate_pipeline_inputs --------------------------------------------

test_that(".validate_pipeline_inputs accepts a well-formed schema and functions", {
  expect_silent(
    .validate_pipeline_inputs(event_schema_se, scorer = identity,
                              opts = splicing_options(),
                              style = splicing_style(), plot_fn = identity)
  )
})

test_that(".validate_pipeline_inputs rejects a schema missing required fields", {
  broken <- event_schema_se
  broken$build_regions <- NULL
  expect_error(
    .validate_pipeline_inputs(broken, identity, splicing_options(),
                              splicing_style(), identity),
    class = "rnapeaks_error_invalid_arg"
  )
})

test_that(".validate_pipeline_inputs rejects non-function scorer / plot_fn", {
  expect_error(
    .validate_pipeline_inputs(event_schema_se, scorer = "nope",
                              splicing_options(), splicing_style(), identity),
    class = "rnapeaks_error_invalid_arg"
  )
  expect_error(
    .validate_pipeline_inputs(event_schema_se, identity, splicing_options(),
                              splicing_style(), plot_fn = 42),
    class = "rnapeaks_error_invalid_arg"
  )
})

# --- .report_group_sizes --------------------------------------------------

test_that(".report_group_sizes prints only in verbose mode", {
  idx <- list(Negative = 1:2, Positive = 1:3, Control = integer(0))
  expect_silent(.report_group_sizes(idx, splicing_options(verbose = FALSE)))
  expect_message(.report_group_sizes(idx, splicing_options(verbose = TRUE)),
                 regexp = "Negative=2")
})

# --- .check_hits_shape ----------------------------------------------------

test_that(".check_hits_shape accepts a valid hit list including the empty case", {
  expect_silent(.check_hits_shape(list(event_id = c(1L, 1L),
                                       col_idx = c(2L, 5L)), n_positions = 6))
  # zero hits is legitimate (a group can score nothing) and skips the range check
  expect_silent(.check_hits_shape(list(event_id = integer(0),
                                       col_idx = integer(0)), n_positions = 6))
})

test_that(".check_hits_shape rejects malformed, mismatched, or out-of-range hits", {
  expect_error(.check_hits_shape(list(event_id = 1L), n_positions = 6),
               class = "rnapeaks_error_invalid_arg")                # missing col_idx
  expect_error(.check_hits_shape(list(event_id = c(1L, 1L),
                                      col_idx = 1L), n_positions = 6),
               class = "rnapeaks_error_invalid_arg")                # length mismatch
  expect_error(.check_hits_shape(list(event_id = 1L, col_idx = 7L),
                                 n_positions = 6),
               class = "rnapeaks_error_invalid_arg")                # col_idx > n
})

# --- .significance_table --------------------------------------------------

test_that(".significance_table returns one significance block per tested group", {
  opts  <- splicing_options(stat_test = "fisher-all", use_fdr = FALSE)
  style <- splicing_style(show_significance = TRUE)
  sig <- .significance_table(per_group_fixture(), opts, style)
  # Negative + Positive tested against Control; Control itself is not tested.
  expect_equal(sort(unique(sig$group)), c("Negative", "Positive"))
  expect_equal(nrow(sig), 2 * 6)                # 2 groups x 6 positions
  expect_true(all(c("position", "pvalue", "significant", "group") %in% names(sig)))
})

test_that(".significance_table skips when disabled, control-less, or control-empty", {
  opts <- splicing_options()
  pg   <- per_group_fixture()
  expect_null(.significance_table(pg, opts, splicing_style(show_significance = FALSE)))
  expect_null(.significance_table(pg[c("Positive", "Negative")], opts,
                                  splicing_style()))
  pg$Control$n <- 0L
  expect_warning(res <- .significance_table(pg, opts, splicing_style()),
                 regexp = "no events")
  expect_null(res)
})

# --- .assemble_frequency_frame --------------------------------------------

test_that(".assemble_frequency_frame divides counts by n and lays out the axis", {
  # n_regions = 2, region_width = 3 -> 6 positions per group.
  pg <- list(
    Positive = list(counts = c(2, 0, 0, 4, 0, 0), n = 4L),
    Control  = list(counts = c(9, 9, 9, 9, 9, 9), n = 3L)   # overwritten below
  )
  cstats <- list(mean_per_position = rep(0.1, 6), sd_per_position = rep(0.2, 6))
  df <- .assemble_frequency_frame(pg, cstats, n_regions = 2, region_width = 3)

  pos <- df[df$group == "Positive", ]
  ctl <- df[df$group == "Control", ]
  expect_equal(pos$frequency, c(0.5, 0, 0, 1, 0, 0))   # counts / 4
  expect_equal(pos$frequency_sd, rep(0, 6))            # only Control carries SD
  expect_equal(ctl$frequency, rep(0.1, 6))             # taken from the bootstrap
  expect_equal(ctl$frequency_sd, rep(0.2, 6))
  # axis geometry: region 1 then region 2, position 1..3 within each
  expect_equal(pos$region_idx, c(1, 1, 1, 2, 2, 2))
  expect_equal(pos$position_in_region, c(1, 2, 3, 1, 2, 3))
  expect_equal(pos$global_position, 1:6)
  expect_equal(unique(pos$n_events), 4L)
})

test_that(".assemble_frequency_frame gives an empty group a flat-zero curve", {
  pg <- list(Positive = list(counts = integer(6), n = 0L))
  df <- .assemble_frequency_frame(pg, NULL, n_regions = 2, region_width = 3)
  expect_equal(df$frequency, rep(0, 6))
  expect_equal(unique(df$n_events), 0L)
})

# --- .smooth_by_group -----------------------------------------------------

test_that(".smooth_by_group smooths each group independently within regions", {
  # One region of width 3 per group. Group A is a spike -> centered width-3 mean
  # c(1.5, 1, 1.5); group B is flat and must stay untouched by A's values.
  values <- c(0, 3, 0, 9, 9, 9)
  groups <- rep(c("A", "B"), each = 3)
  out <- .smooth_by_group(values, groups, window = 3,
                          n_regions = 1, region_width = 3)
  expect_equal(out, c(1.5, 1, 1.5, 9, 9, 9))
})

# --- prepare_event_map ----------------------------------------------------

se_events <- function() {
  data.frame(
    chr = "chr1", strand = "+", GeneID = c("G1", "G2", "G3"),
    PValue = c(0.01, 0.01, 0.99),
    FDR    = c(0.01, 0.01, 0.99),
    IncLevelDifference = c(0.30, -0.30, 0.001),   # Positive, Negative, Control
    exonStart_0base = 1400, exonEnd = 1500,
    upstreamES = 1000, upstreamEE = 1100,
    downstreamES = 1800, downstreamEE = 1900,
    stringsAsFactors = FALSE
  )
}

test_that("prepare_event_map normalizes chr, partitions groups, and builds regions", {
  opts <- splicing_options(min_count = 0, width_exon = 10, width_intron = 20)
  prep <- suppressMessages(prepare_event_map(se_events(), event_schema_se, opts))

  expect_equal(prep$n_regions, 4L)
  expect_equal(prep$region_width, 31L)             # 10 + 20 + 1
  expect_equal(prep$total_positions, 124L)         # 4 * 31
  expect_equal(unique(prep$events$chr), "1")       # "chr1" -> "1"
  expect_equal(prep$groups_idx$Positive, 1L)
  expect_equal(prep$groups_idx$Negative, 2L)
  expect_equal(prep$groups_idx$Control, 3L)
  # one event in each group -> 4 region ranges apiece
  expect_length(prep$regions_by_group$Positive, 4L)
  expect_s4_class(prep$regions_by_group$Positive, "GRanges")
  expect_null(prep$schematic_data)                 # SE has no schematic_stat
})

test_that("prepare_event_map computes schematic stats when the schema defines them", {
  opts <- splicing_options(min_count = 0)
  # Two Positive events with long-isoform extensions of 100 and 200 bp; the
  # empty Negative / Control groups warn (tested in test-sm-filter) and are
  # irrelevant to the schematic stat, so those warnings are suppressed here.
  a5 <- data.frame(
    chr = "1", strand = "+", GeneID = c("G1", "G2"),
    PValue = 0.01, FDR = 0.01, IncLevelDifference = c(0.3, 0.4),
    longExonStart_0base = 1000, longExonEnd = c(1200, 1300),   # long-short = 100, 200
    shortES = 1000, shortEE = 1100,
    flankingES = 1500, flankingEE = 1600,
    stringsAsFactors = FALSE
  )
  prep <- suppressWarnings(suppressMessages(
    prepare_event_map(a5, event_schema_a5ss, opts)
  ))
  expect_equal(prep$schematic_data$median_ext, 150)   # median(c(100, 200))
})

# --- event_map_pipeline ---------------------------------------------------

test_that("event_map_pipeline scores, tabulates, and hands a frequency frame to plot_fn", {
  # Deterministic stub scorer: fixed hits per group so counts are predictable.
  scorer <- function(regions, n_regions, region_width, group_name) {
    if (identical(group_name, "Positive")) {
      list(event_id = c(1L, 1L), col_idx = c(1L, 4L))   # counts hit pos 1 and 4
    } else {
      list(event_id = 1L, col_idx = 1L)                 # Control hits pos 1
    }
  }
  captured <- new.env()
  plot_fn <- function(data, schema, style, opts, significance, title,
                      schematic_data) {
    captured$data <- data
    captured$sig  <- significance
    "PLOT"
  }

  prep <- list(
    groups_idx       = list(Positive = 1:5, Control = 1:3),
    events           = data.frame(id = 1:5),
    regions_by_group = list(Positive = GenomicRanges::GRanges(),
                            Control  = GenomicRanges::GRanges()),
    n_regions        = 2L,
    region_width     = 3L,
    total_positions  = 6L,
    schematic_data   = NULL
  )
  opts <- splicing_options(stat_test = "fisher-all", moving_average = 0)

  set.seed(1)
  res <- suppressMessages(
    event_map_pipeline(schema = event_schema_ri, scorer = scorer, opts = opts,
                       style = splicing_style(), plot_fn = plot_fn, prep = prep)
  )

  expect_equal(res$plot, "PLOT")
  # Positive frequency = counts / n = c(1,0,0,1,0,0) / 5
  pos <- captured$data[captured$data$group == "Positive", ]
  expect_equal(pos$frequency, c(0.2, 0, 0, 0.2, 0, 0))
  expect_true("moving_avg" %in% names(captured$data))
  # significance computed for Positive vs Control across all 6 positions
  expect_equal(unique(captured$sig$group), "Positive")
  expect_equal(nrow(captured$sig), 6)
  # the Positive events are returned for downstream use; Negative was not scored
  expect_equal(nrow(res$data$Positive), 5)
  expect_null(res$data$Negative)
})

test_that("event_map_pipeline requires either events or prep", {
  expect_error(
    event_map_pipeline(schema = event_schema_se, scorer = identity,
                       opts = splicing_options(), style = splicing_style(),
                       plot_fn = identity),
    class = "rnapeaks_error_invalid_arg"
  )
})
