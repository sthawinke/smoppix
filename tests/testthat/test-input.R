context("Test input smoppix package")
test_that("Reading in data proceeds without errors", {
    expect_message(buildHyperFrame(df, coordVars = c("x", "y"), imageVars = "fov"))
    # Should work also with single image variable
    expect_message(hypFrame <- buildHyperFrame(df, coordVars = c("x", "y"), imageVars = c(
        "fov", "condition", "age"
    ), imageIdentifier = c("fov", "condition")))
    expect_s3_class(hypFrame, c("hyperframe", "list"))
    expect_message(buildHyperFrame(as.matrix(df[, c("x", "y")]),
        imageVars = df$fov,
        covariates = df[, c("gene", "condition", "age"), drop = FALSE]
    ))
    expect_message(hypFrame2 <- buildHyperFrame(as.matrix(df[, c("x", "y")]), imageIdentifier = df[,c("fov", "condition")], imageVars = df[
        ,c("fov", "condition", "age")
    ], covariates = df[, "gene", drop = FALSE]))
    expect_identical(hypFrame, hypFrame2)
    expect_silent(hypFrame3 <- buildHyperFrame(lapply(listPPP, identity)))
    expect_silent(hypFrame4 <- buildHyperFrame(split(df, f = apply(df[, c("fov","condition")], 1, paste, collapse = "_")), 
                                                covariatesDf = df[!duplicated(df[,c("fov", "condition")]), c("fov", "condition")]))
})
# Read in spatial experiment
library(SpatialExperiment)
library(DropletUtils)
example(read10xVisium)
test_that("Reading in SpatialExperiment class proceeds without errors", {
    expect_message(hypFrame4 <- buildHyperFrame(spe, imageVars = "sample_id", pointVars = "in_tissue"))
})
test_that("Adding regions of interest works", {
    expect_s3_class(hypFrame5 <- addCell(hypFrame, wList), "hyperframe")
    expect_s3_class(
        hypFrame6 <- addCell(hypFrame, wList, cellTypes = cellTypesDf),
        "hyperframe"
    )
    expect_error(addCell(hypFrame, wList2, findOverlappingOwins = TRUE)) # Detect overlap
    expect_type(marks(hypFrame5[[1, "ppp"]])$cell, "character")
    expect_type(marks(hypFrame6[[1, "ppp"]])$cellType, "character")
    expect_error(hypFrame4 <- buildHyperFrame(spe,
        imageVars = c("sample_id", "in_tissue"),
        coVars = "in_tissue"
    ))
})

test_that("Adding regions of interest throws errors when appropriate", {
    wListNoNames <- wList
    names(wListNoNames) <- NULL
    expect_error(addCell(addCell(hypFrame, wList), wList))
    expect_error(addCell(hypFrame$ppp[[1]], wList))
    expect_error(addCell(hypFrame, wListNoNames))
})
