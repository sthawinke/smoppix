context("Constructing the weight function")
test_that("Adding the weight function works", {
    expect_false(is.null((objBG <- addWeightFunction(piEstsBG, designVars = "condition"))$Wfs))
    expect_false(is.null((objCSR <- addWeightFunction(piEstsCSR, designVars = "condition"))$Wfs))
    # Alternatively, specifying the lowest level
    expect_false(is.null((objBG2 <- addWeightFunction(piEstsBG, lowestLevelVar = "fov"))$Wfs))
    expect_silent(objCSR2 <- addWeightFunction(piEstsCSR, pi = c("nnCell", "nnCellPair")))
    expect_true(all(vapply(objBG$wfs, FUN.VALUE = TRUE, is, "scam")))
    expect_true(all(vapply(objCSR$wfs, FUN.VALUE = TRUE, is, "scam")))
    expect_type(pred <- evalWeightFunction(objBG$Wfs[["nn"]]), "double")
    expect_true(all(pred > 0))
    expect_type(
        predNew <- evalWeightFunction(objBG$Wfs[["nn"]], newdata = data.frame(NP = 5)),
        "double"
    )
    expect_true(predNew > evalWeightFunction(objBG$Wfs[["nn"]], newdata = data.frame(NP = 2)))
    # With information sharing across features
    expect_s3_class(dfBG1 <- buildDataFrame(objBG, gene = "gene1", pi = "nn"), "data.frame")
    expect_s3_class(
        dfBG2 <- buildDataFrame(objBG, gene = "gene1--gene2", pi = "nnPair"),
        "data.frame"
    )
    expect_warning(addWeightFunction(objBG2, lowestLevelVar = "fov"))
})
test_that("Weight function application throws errors where appropriate", {
    expect_error(buildDataFrame(piEstsBG, gene = "gene101", pi = "nn"))
    expect_error(buildDataFrame(piEstsBG, gene = "gene1--gene2", pi = "nnPair"))
    expect_error(buildDataFrame(piEstsBG, gene = "gene1", pi = "nn"))
    expect_error(addWeightFunction(piEstsCSR, designVars = "condition", lowestLevelVar = "fov"))
    expect_error(addWeightFunction(piEstsCSR, pi = "nn", lowestLevelVar = c(
        "fov",
        "condition"
    )))
    expect_error(addWeightFunction(piEstsCSR, pi = "nn", lowestLevelVar = "treatment"))
    expect_error(addWeightFunction(hypYang))
})
