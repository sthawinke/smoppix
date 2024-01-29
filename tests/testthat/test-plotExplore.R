context("Test plotExplore function")
test_that("Plotting hyperframes proceeds without errors", {
    expect_silent(plotExplore(hypYang))
    expect_silent(plotExplore(hypFrame2))
    expect_silent(plotExplore(hypFrame2, piEsts = objCSR, piColourCell = "edge", feature = "gene1"))
})

test_that("Plotting hyperframes throws warning when appropriate", {
    expect_warning(plotExplore(hypYang, plotWindows = TRUE))
    expect_error(plotExplore(hypFrame2, piEsts = objCSR, piColourCell = "edge"))
})
