# Tests for filter_events() / .apply_count_filter() in R/sm-event_filter.R
#
# Group partitioning rules (trio schemas):
#   sig      = FDR < event_fdr  AND PValue < event_fdr
#   control  = PValue > control_pval AND FDR > control_pval
#   Negative = sig & IncLevelDifference < psi_cutoff[1]
#   Positive = sig & IncLevelDifference > psi_cutoff[2]
#   Control  = control & abs(IncLevelDifference) < psi_control_max
# Defaults: event_fdr .05, control_pval .95, psi_cutoff c(-.1, .1),
#           psi_control_max .005.

# One event engineered for each outcome.
partition_events <- function() {
  data.frame(
    FDR                = c(0.01, 0.01, 0.99, 0.50, 0.01),
    PValue             = c(0.01, 0.01, 0.99, 0.50, 0.01),
    IncLevelDifference = c(0.30, -0.30, 0.001, 0.02, 0.05),
    stringsAsFactors = FALSE
  )
  # row1 Positive, row2 Negative, row3 Control, row4 none, row5 sig but |dPSI|<cutoff -> none
}

test_that("events are partitioned into Negative / Positive / Control by the rules", {
  opts <- splicing_options(min_count = 0)   # count filter off
  res <- filter_events(partition_events(), event_schema_se, opts)
  expect_equal(res$groups_idx$Positive, 1L)
  expect_equal(res$groups_idx$Negative, 2L)
  expect_equal(res$groups_idx$Control,  3L)
})

test_that("only the requested groups are returned, in canonical order", {
  opts <- splicing_options(min_count = 0, groups = c("Positive", "Control"))
  res <- filter_events(partition_events(), event_schema_se, opts)
  expect_equal(names(res$groups_idx), c("Positive", "Control"))
})

test_that("a requested group with no events warns", {
  opts <- splicing_options(min_count = 0)
  no_control <- partition_events()[1:2, ]   # Positive + Negative only
  expect_warning(filter_events(no_control, event_schema_se, opts),
                 regexp = "no events")
})

test_that("a Single-group schema returns every row in one group", {
  opts <- splicing_options(min_count = 0)
  res <- filter_events(partition_events(), list(group_set = "Single"), opts)
  expect_equal(names(res$groups_idx), "Single")
  expect_equal(res$groups_idx$Single, seq_len(nrow(partition_events())))
})

# --- .apply_count_filter --------------------------------------------------

count_events <- function() {
  data.frame(
    IJC_SAMPLE_1 = c("5,5", "1"),   # summed replicates: 10, 1
    SJC_SAMPLE_1 = c("10",  "2"),   #                    10, 2
    IJC_SAMPLE_2 = c("8",   "50"),
    SJC_SAMPLE_2 = c("8",   "50"),
    stringsAsFactors = FALSE
  )
}

test_that("count filter keeps events whose summed junction counts exceed min_count", {
  opts <- splicing_options(min_count = 10)
  out <- .apply_count_filter(count_events(), opts)
  # row1: sample1 = 20, sample2 = 16 (both > 10) -> kept
  # row2: sample1 = 3 (<= 10) -> dropped
  expect_equal(nrow(out), 1)
})

test_that("min_count = 0 disables the count filter entirely", {
  opts <- splicing_options(min_count = 0)
  expect_equal(nrow(.apply_count_filter(count_events(), opts)), 2)
})

test_that("count filter aborts when the junction-count columns are missing", {
  opts <- splicing_options(min_count = 10)
  expect_error(.apply_count_filter(data.frame(x = 1:2), opts),
               class = "rnapeaks_error_invalid_arg")
})
