context("Constructing the weight function")

test_that("Adding the weight function works", {
    expect_silent(objBG <- addWeightFunction(piEstsBG, designVars = "condition"))
    expect_silent(objCSR <- addWeightFunction(piEstsCSR, designVars = "condition"))
    #Alternatively, specifying the lowest level
    expect_silent(objBG2 <- addWeightFunction(piEstsBG, lowestLevelVar = "fov"))
    expect_silent(objCSR2 <- addWeightFunction(piEstsCSR, lowestLevelVar = "fov"))
    expect_equal(objBG, objBG2, tolerance = 0.4)
    #Tolerate differences in CPU time
    expect_equal(objCSR, objCSR2, tolerance = 0.4)
    expect_true(all(vapply(objBG$wfs, FUN.VALUE = TRUE, is, "scam")))
    expect_true(all(vapply(objCSR$wfs, FUN.VALUE = TRUE, is, "scam")))
    expect_type(pred <- evalWeightFunction(objBG$Wfs[["nn"]]), "double")
    expect_true(all(pred>0))
    expect_type(predNew <- evalWeightFunction(objBG$Wfs[["nn"]], newdata = data.frame("NP" = 5)), "double")
    expect_true(predNew > evalWeightFunction(objBG$Wfs[["nn"]], newdata = data.frame("NP" = 4)))
    #With information sharing across features
    expect_s3_class(dfBG1 <- buildDfMM(objBG, gene = "gene1", pi = "nn"), "data.frame")
    expect_s3_class(dfBG2 <- buildDfMM(objBG, gene = "gene1--gene2", pi = "nnPair"),
                    "data.frame")
})
test_that("Weight function application throws errors where appropriate", {
    expect_error(buildDfMM(piEstsBG, gene = "gene101", pi = "nn"))
    expect_error(buildDfMM(piEstsBG, gene = "gene1--gene2", pi = "nnPair"))
    expect_error(buildDfMM(piEstsBG, gene = "gene1", pi = "nn"))
    expect_error(addWeightFunction(piEstsCSR, designVars = "condition", lowestLevelVar = "fov"))
})
