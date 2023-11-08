context("Constructing the weight function")
test_that("The built weight function has the desired properties", {
    expect_s3_class(wf <- buildWeightFunction(piEstsBG, pi = "nn", hypFrame = hypFrame2, designVars = "condition"), "scam")
    expect_s3_class(wfCSR <- buildWeightFunction(piEstsCSR, pi = "nn", hypFrame = hypFrame2, designVars = "condition"), "scam")
    expect_type(pred <- evalWeightFunction(wf), "numeric")
    expect_true(all(pred>0))
    expect_type(predNew <- evalWeightFunction(wf, newdata = data.frame("NP" = 5)), "numeric")
    expect_true(predNew > evalWeightFunction(wf, newdata = data.frame("NP" = 4)), "numeric")
})
