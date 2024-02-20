context("Writing results to spreadsheets")
test_that("Sheet naming works",{
    expect_identical(makeSheetName("edge", "Intercept", TRUE), "Close to edge | overall")
    expect_identical(makeSheetName("centroid", "Intercept", FALSE), "Far from centroid | overall")
    expect_identical(makeSheetName("nn", "day"), "Aggregated | day")
    expect_identical(makeSheetName("nnPair", "Intercept", FALSE), "Antilocalized | overall")
})
