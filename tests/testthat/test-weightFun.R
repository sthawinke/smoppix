context("Constructing the weight function")
test_that("The built weight function has the desired properties", {
    expect_s3_class(wf <- buildWeightFunction(piEstsBG, pi = "nn", hypFrame = hypFrame2, designVars = "condition"), "scam")
})
