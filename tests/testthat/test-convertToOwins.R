context("Test conversion to owins")
test_that("Conversion of rois works", {
})
test_that("Conversion of SpatialPolygonsDataFrame works", {
})
# Regions of interest (roi): Diamond in the center plus four triangles
w1 <- data.frame(x = c(0, .5, 1, .5), y = c(.5, 0, .5, 1))
w2 <- data.frame(x = c(0, 0, .5), y = c(.5, 0, 0))
w3 <- data.frame(x = c(0, 0, .5), y = c(1, 0.5, 1))
w4 <- data.frame(x = c(1, 1, .5), y = c(0.5, 1, 1))
w5 <- data.frame(x = c(1, 1, .5), y = c(0, 0.5, 0))
wListDf <- lapply(seq_len(nDesignFactors), function(x) {
    Li <- list("w1" = w1, "w2" = w2, "w3" = w3, "w4" = w4, "w5" = w5)
    names(Li) <- paste0(names(Li), "_", x)
    Li
})
names(wListDf) <- rownames(hypFrame)
test_that("Conversion of coordinate matrices works", {
    expect_s3_class(acDf <- addCell(hypFrame, wListDf), "hyperframe")
    expect_is(marks(acDf$ppp[[1]])$cell, "character")
})
