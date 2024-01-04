context("Adding and plotting the weight function")
test_that("Plotting the weight function works", {
    expect_silent(plotWf(yangPims, "nn"))
    expect_silent(plotWf(yangPims, "nnPair"))
    expect_silent(plotWf(objBG, "nn"))
    expect_silent(plotWf(objBG, "nnCell"))
    expect_silent(plotWf(objBG, "nnPair"))
    expect_silent(plotWf(objBG, "nnPairCell"))
})
test_that("Plotting the weight function application throws errors where appropriate", {
    expect_error(plotWf(piEstsBG))
})
