context("Constructing the weight function")
test_that("The built weight function has the desired properties", {
    expect_true(all(vapply(objBG$wfs, FUN.VALUE = TRUE, is, "scam")))
    expect_true(all(vapply(objCSR$wfs, FUN.VALUE = TRUE, is, "scam")))
    expect_s3_class(wfPairAll <- buildWeightFunction(piEstsBG, pi = "allDistPair", hypFrame = hypFrame2, designVars = "condition"), "scam")
    expect_type(pred <- evalWeightFunction(objBG$wfs[["nn"]]), "double")
    expect_true(all(pred>0))
    expect_type(predNew <- evalWeightFunction(wf, newdata = data.frame("NP" = 5)), "double")
    expect_true(predNew > evalWeightFunction(wf, newdata = data.frame("NP" = 4)))
    #With information sharing across features
    expect_s3_class(dfBG1 <- buildDfMM(piEstsBG, gene = "gene1", pi = "nn", hypFrame = hypFrame2, weightFunction = wf), "data.frame")
    expect_s3_class(dfBG2 <- buildDfMM(piEstsBG, gene = "gene1--gene2", pi = "nnPair", hypFrame = hypFrame2, weightFunction = wfPair),
                    "data.frame")
})
test_that("Weight function application throws errors where appropriate", {
    wf <- buildWeightFunction(piEstsBG, pi = "nn", hypFrame = hypFrame2, designVars = "condition")
    expect_error(buildDfMM(piEstsBG, gene = "gene101", pi = "nn", hypFrame = hypFrame2, weightFunction = wf))
    expect_error(buildDfMM(piEstsBG, gene = "gene1--gene2", pi = "nnPair", hypFrame = hypFrame2, weightFunction = wf))
    expect_error(buildDfMM(piEstsBG, gene = "gene1", pi = "nn", hypFrame = hypFrame2, weightFunction = wfPair))
})
