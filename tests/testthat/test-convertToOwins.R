context("Test conversion to owins")
test_that("Conversion of rois works", {
    expect_message(buildHyperFrame(df, coordVars = c("x", "y"),
                                   imageVars = "fov"))
    # Should work also with single image variable

})
test_that("Conversion of SpatialPolygonsDataFrame works", {
    expect_s3_class(hypFrame5 <- addCell(hypFrame, wList), "hyperframe")
    expect_s3_class(hypFrame6 <- addCell(hypFrame, wList,
                                         cellTypes = cellTypesDf), "hyperframe")
    expect_warning(addCell(hypFrame, wList2)) # Detect overlap
    expect_type(marks(hypFrame5[[1, "ppp"]])$cell, "character")
    expect_type(marks(hypFrame6[[1, "ppp"]])$cellType, "character")
})
# Regions of interest (roi): Diamond in the center plus four triangles
w1 <- data.frame(x = c(0, .5, 1, .5), y = c(.5, 0, .5, 1))
w2 <- data.frame(x = c(0, 0, .5), y = c(.5, 0, 0))
w3 <- data.frame(x = c(0, 0, .5), y = c(1, 0.5, 1))
w4 <- data.frame(x = c(1, 1, .5), y = c(0.5, 1, 1))
w5 <- data.frame(x = c(1, 1, .5), y = c(0, 0.5, 0))
test_that("Conversion of coordinate matrices works", {
    expect_s3_class(addCell(hypFrame, wList), "hyperframe")
    expect_error(addCell(hypFrame$ppp[[1]], wList))
})
