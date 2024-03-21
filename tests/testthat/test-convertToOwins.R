context("Test conversion to owins")
if (requireNamespace("RImageJROI")) {
    path <- file.path(system.file(package = "RImageJROI"), "extdata", "ijroi")
    rect <- RImageJROI::read.ijroi(file.path(path, "rect.roi"))
    poly <- RImageJROI::read.ijroi(file.path(path, "polygon.roi"))
    oval <- RImageJROI::read.ijroi(file.path(path, "oval.roi"))
    wListRoi <- lapply(seq_len(nDesignFactors), function(x) {
        Li <- list(w1 = rect, w2 = poly, w3 = oval)
        names(Li) <- paste0(names(Li), "_", x)
        Li
    })
    names(wListRoi) <- rownames(hypFrame)
    hypFrameRoi <- hypFrame
    for (x in seq_len(nrow(hypFrameRoi))) {
        PPP <- hypFrameRoi$ppp[[x]]
        PPP$window <- owin(xrange = c(0, 100), yrange = c(0, 100))
        coords(PPP) <- data.frame(x = runif(npoints(PPP), 1, 100), y = runif(
            npoints(PPP),
            1, 100
        ))
        hypFrameRoi$ppp[[x]] <- PPP
    }
    test_that("Conversion of rois works", {
        expect_s3_class(
            acRoi <- addCell(hypFrameRoi, wListRoi, warnOut = FALSE),
            "hyperframe"
        )
        expect_is(marks(acRoi$ppp[[1]])$cell, "character")
    })
}
if (requireNamespace("polyCub")) {
    library(sp)
    grd <- GridTopology(c(1, 1), c(1, 1), c(10, 10))
    polys <- as(grd, "SpatialPolygons")
    centroids <- coordinates(polys)
    x <- centroids[, 1]
    y <- centroids[, 2]
    z <- 1.4 + 0.1 * x + 0.2 * y + 0.002 * x * x
    ex_1.7 <- SpatialPolygonsDataFrame(polys, data = data.frame(
        x = x, y = y, z = z,
        row.names = row.names(polys)
    ))
    wListSPDF <- lapply(seq_len(nDesignFactors), function(x) {
        ex_1.7
    })
    names(wListSPDF) <- rownames(hypFrame)
    test_that("Conversion of SpatialPolygonsDataFrame works", {
        expect_s3_class(
            acSPDF <- addCell(hypFrame, wListSPDF, warnOut = FALSE),
            "hyperframe"
        )
        expect_is(marks(acSPDF$ppp[[1]])$cell, "character")
    })
}
if (requireNamespace("polyCub")) {
    # Regions of interest (roi): Diamond in the center plus four triangles
    w1 <- data.frame(x = c(0, 0.5, 1, 0.5), y = c(0.5, 0, 0.5, 1))
    w2 <- data.frame(x = c(0, 0, 0.5), y = c(0.5, 0, 0))
    w3 <- data.frame(x = c(0, 0, 0.5), y = c(1, 0.5, 1))
    w4 <- data.frame(x = c(1, 1, 0.5), y = c(0.5, 1, 1))
    w5 <- data.frame(x = c(1, 1, 0.5), y = c(0, 0.5, 0))
    wListDf <- lapply(seq_len(nDesignFactors), function(x) {
        Li <- list(w1 = w1, w2 = w2, w3 = w3, w4 = w4, w5 = w5)
        names(Li) <- paste0(names(Li), "_", x)
        Li
    })
    names(wListDf) <- rownames(hypFrame)
    test_that("Conversion of coordinate matrices works", {
        expect_s3_class(acDf <- addCell(hypFrame, wListDf), "hyperframe")
        expect_is(marks(acDf$ppp[[1]])$cell, "character")
    })
}
test_that("Error is thrown in case of wrong window type", {
    expect_error(addCell(hypFrame, list(w1 = rnorm, w2 = rpois)))
})
