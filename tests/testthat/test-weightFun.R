context("Constructing the weight function")
test_that("The built weight function has the desired properties", {
    expect_s3_class(wf <- buildWeightFunction(piEstsBG, pi = "nn", hypFrame = hypFrame2, designVars = "condition"), "scam")
    expect_s3_class(wfCSR <- buildWeightFunction(piEstsCSR, pi = "nn", hypFrame = hypFrame2, designVars = "condition"), "scam")
    expect_s3_class(wfPair <- buildWeightFunction(piEstsBG, pi = "nnPair", hypFrame = hypFrame2, designVars = "condition"), "scam")
    expect_s3_class(wfAll <- buildWeightFunction(piEstsBG, pi = "allDist", hypFrame = hypFrame2, designVars = "condition"), "scam")
    expect_s3_class(wfPairAll <- buildWeightFunction(piEstsBG, pi = "allDistPair", hypFrame = hypFrame2, designVars = "condition"), "scam")
    expect_type(pred <- evalWeightFunction(wf), "double")
    expect_true(all(pred>0))
    expect_type(predNew <- evalWeightFunction(wf, newdata = data.frame("NP" = 5)), "double")
    expect_true(predNew > evalWeightFunction(wf, newdata = data.frame("NP" = 4)))
})
test_that("The weight function can be used to build the desired dataframe", {
    expect_s3_class(dfBG1 <- buildDfMM(piEstsBG, gene = "gene1", pi = "nn", hypFrame = hypFrame2, weightFunction = wf), "data.frame")
    expect_s3_class(dfBG2 <- buildDfMM(piEstsBG, gene = "gene1_gene2", pi = "nnPair", hypFrame = hypFrame2, weightFunction = wfPair),
                    "data.frame")
})
