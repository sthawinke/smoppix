context("Test plotCell function")
test_that("Plotting cells proceeds without errors", {
    expect_silent(plotCells(hypFrame2, "gene1"))
    expect_silent(plotCells(hypFrame2, "gene1", borderColVar = "condition"))
    expect_silent(plotCells(hypFrame2, "gene1", borderColVar = "cellType", 
                            scaleBarSize = c(0.01, 0.1)))
})

test_that("Plotting cells throws error when appropriate", {
    expect_error(plotCells(hypFrame2, "gene1", borderColVar = "bar"))
    expect_error(plotCells(hypFrame2, "gene1", scaleBarSize = c(1,20, 9)))
})
