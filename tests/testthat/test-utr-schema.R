# Tests for event_schema_utr in R/utr-schema.R

test_that("the schema declares a 100-bin window", {
  expect_equal(event_schema_utr$n_bins, 100L)
})

test_that("plot_layout returns the region / bin / cds geometry", {
  lay <- event_schema_utr$plot_layout()
  expect_equal(lay$region_start, 0L)
  expect_equal(lay$bin_width, 100L)
  expect_equal(lay$cds_width, 25L)

  custom <- event_schema_utr$plot_layout(n_bins = 50L, cds_width = 10L)
  expect_equal(custom$bin_width, 50L)
  expect_equal(custom$cds_width, 10L)
})

test_that("build_schematic_layers puts the CDS block on the side-appropriate end", {
  lay <- event_schema_utr$plot_layout()
  layers5 <- event_schema_utr$build_schematic_layers(lay, utr_style(),
                                                     y_min = 0, exon_height = 1,
                                                     side = "utr5")
  layers3 <- event_schema_utr$build_schematic_layers(lay, utr_style(),
                                                     y_min = 0, exon_height = 1,
                                                     side = "utr3")
  expect_length(layers5, 5)
  expect_s3_class(layers5[[1]], "LayerInstance")
  # CDS box (first layer) sits to the RIGHT of a 5' UTR and LEFT of a 3' UTR
  expect_equal(layers5[[1]]$data$xmin, 100)   # bin_width
  expect_equal(layers3[[1]]$data$xmin, -25)   # -cds_width
})
