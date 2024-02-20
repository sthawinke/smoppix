context("Writing results to spreadsheets")
test_that("Sheet naming works",{
    expect_identical(makeSheetName("edge", "Intercept", TRUE), "Close to edge|baseline")
    expect_identical(makeSheetName("centroid", "Intercept", FALSE), "Far from centroid|baseline")
    expect_identical(makeSheetName("nn", "day"), "Aggregation|day")
    expect_identical(makeSheetName("nnCell", "veryLongVariable"), "Aggregation in cell|veryLongVar")
    expect_identical(makeSheetName("nnPair", "Intercept", FALSE), "Antilocalization|baseline")
})
