context("Plotting the weight function")
wf <- buildWeightFunction(piEstsBG, pi = "nn", hypFrame = hypFrame2, designVars = "condition")
wfPair <- buildWeightFunction(piEstsBG, pi = "nnPair", hypFrame = hypFrame2, designVars = "condition")
test_that("Plotting the weight function works", {
    expect_s3_class(plotWf(wf), "plot")
    expect_s3_class(plotWf(wf), "ggplot")
})
test_that("Plotting the weight function application throws errors where appropriate", {
    expect_error(plotWf(piEstsBG))
})
