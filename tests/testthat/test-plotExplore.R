context("Test plotExplore function")
test_that("Plotting hyperframes proceeds without errors", {
    expect_silent(plotExplore(hypYang))
    expect_silent(plotExplore(hypFrame2))
    expect_silent(plotExplore(hypFrame2, piEsts = objCSR, piColourCell = "edge",
        feature = "gene1"))
    expect_silent(plotExplore(hypFrame2, piEsts = objCSR, piColourCell = "nnCell",
        feature = "gene1"))
    expect_silent(plotExplore(hypFrame2, piEsts = objCSR, piColourCell = "nnPairCell",
        feature = "gene1--gene2"))
})

test_that("Plotting hyperframes throws warning when appropriate", {
    expect_warning(plotExplore(hypYang, plotWindows = TRUE))
    expect_error(plotExplore(hypFrame2, piEsts = objCSR, piColourCell = "edge"))
    expect_error(plotExplore(hypFrame2, piEsts = objCSR, piColourCell = "nn", feature = "gene1"))
    expect_error(plotExplore(hypFrame2, piEsts = objCSR, piColourCell = "nnPairCell",
        feature = "gene2"))
    expect_error(plotExplore(hypFrame2, piEsts = objCSR, piColourCell = "nnCell",
        feature = "gene1--gene2"))
})
