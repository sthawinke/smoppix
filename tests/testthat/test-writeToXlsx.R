context("Writing results to spreadsheets")
test_that("Sheet naming works", {
    expect_identical(makeSheetName("edge", "Intercept", TRUE), "Close to edge_baseline")
    expect_identical(makeSheetName("centroid", "Intercept", FALSE), "Far from centroid_baseline")
    expect_identical(makeSheetName("nn", "day"), "Aggregation_day")
    expect_identical(makeSheetName("nnCell", "veryLongVariable"), "Aggregation in cell_veryLongVar")
    expect_identical(makeSheetName("nnPair", "Intercept", FALSE), "Antilocalized_baseline")
})
test_that("Writing to xlsx works", {
    tmpFile <- "tmpFile.xlsx"
    fits <- fitLMMs(yangPims,
        pi = "nn", verbose = FALSE,
        features = getFeatures(yangPims)[1:5]
    )
    writeToXlsx(fits, file = tmpFile)
    expect_true(file.exists(tmpFile))
    file.remove(tmpFile)
})
