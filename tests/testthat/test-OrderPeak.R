test_that("OrderPeak orders by count correctly", {
  bed <- data.frame(
    group_name = c("RBFOX2", "RBFOX2", "RBFOX2", "hnRNPC", "hnRNPC", "U2AF2"),
    start = c(100, 200, 300, 100, 200, 100),
    end = c(150, 250, 350, 150, 250, 150),
    stringsAsFactors = FALSE
  )

  result <- OrderPeak(bed, order_by = "Count")

  # RBFOX2 has 3 peaks, hnRNPC has 2, U2AF2 has 1
  expect_equal(result[1], "RBFOX2")
  expect_equal(result[2], "hnRNPC")
  expect_equal(result[3], "U2AF2")
})

test_that("OrderPeak orders by name correctly", {
  bed <- data.frame(
    group_name = c("RBFOX2", "hnRNPC", "U2AF2"),
    start = c(100, 100, 100),
    end = c(150, 150, 150),
    stringsAsFactors = FALSE
  )

  result <- OrderPeak(bed, order_by = "Target")

  # Alphabetical order
  expect_equal(result[1], "RBFOX2")
  expect_equal(result[2], "U2AF2")
  expect_equal(result[3], "hnRNPC")
})

test_that("OrderPeak errors on invalid order_by",
 {
  bed <- data.frame(
    group_name = c("RBFOX2"),
    start = c(100),
    end = c(150),
    stringsAsFactors = FALSE
  )

  expect_error(
    OrderPeak(bed, order_by = "Invalid"),
    "Please provide acceptable feature to rank by"
  )
})
