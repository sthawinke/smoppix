context("Adding and plotting the weight function")
test_that("Plotting the weight function works", {
    expect_silent(plotWf(yangPims, "nn"))
    expect_is(plotWf(yangPims, "nnPair"), "ggplot")
    expect_silent(plotWf(objBG, "nn"))
    expect_silent(plotWf(objBG, "nnCell"))
    expect_is(plotWf(objBG, "nnPair"), "ggplot")
    expect_is(plotWf(objBG, "nnPairCell"), "ggplot")
})
test_that("Plotting the weight function application throws errors where appropriate", {
    expect_error(plotWf(piEstsBG))
})
