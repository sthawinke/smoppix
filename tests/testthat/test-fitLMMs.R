context("Large scale mixed model fitting")
test_that("Fitting linear mixed models proceeds without errors", {
    expect_is(linModsNNint <- fitLMMs(yangPims, pi = "nn", features = getFeatures(yangPims)[1:5]), "list")
    expect_is(fitLMMs(yangPims, pi = "nn", verbose = FALSE, features = getFeatures(yangPims)[1:5]), "list")
    # Supply your own formula
    expect_is(fitLMMs(yangPims, features = getFeatures(yangPims)[1:5],
        pi = "nn",
        Formula = "pi - 0.5 ~ day +1|root"
    ), "list")
    expect_is(linModsNNPairint <- fitLMMs(yangPims,
        features = getFeatures(yangPims)[1:5], pi = "nnPair"), "list")
    expect_is(linModsNN <- fitLMMs(yangPims,
        features = getFeatures(yangPims)[1:5], fixedVars = "day", pi = "nn"),
        "list")
    expect_is(linModsNNPair <- fitLMMs(yangPims,
        features = getFeatures(yangPims)[1:5],
        fixedVars = "day", pi = "nnPair"), "list")
    expect_is(linMModsNN <- fitLMMs(yangPims, features = getFeatures(yangPims)[1:5],
        fixedVars = "day",
        randomVars = "root", pi = "nn"
    ), "list")
    expect_is(linMModsNNPair <- fitLMMs(yangPims,
        fixedVars = "day", features = getFeatures(yangPims)[1:5],
        randomVars = "root", pi = "nnPair"
    ), "list")
    # Returning the models
    expect_is(linModsNNfull <- fitLMMs(yangPims, features = getFeatures(yangPims)[1:5],
        fixedVars = "day",
        pi = "nn", returnModels = TRUE
    ), "list")
    expect_is(linMModsNNfull <- fitLMMs(yangPims, features = getFeatures(yangPims)[1:5],
        fixedVars = "day",
        randomVars = "root", pi = "nn", returnModels = TRUE
    ), "list")
    expect_s3_class(linModsNNfull$models[[1]], "lm")
    expect_s4_class(linMModsNNfull$models[[1]], "lmerModLmerTest")
    expect_is(linModsMP <- fitLMMs(objBG, features = getFeatures(objBG)[1:5],
        fixedVars = "condition",
        pi = "midpoint"
    ), "list")
    expect_is(linModsEdge <- fitLMMs(objBG, features = getFeatures(objBG)[1:5],
        fixedVars = "condition",
        pi = "edge"
    ), "list")
    # Including cell (type) either as fixed or random effect
    # E.g. test for differences between cell types
    expect_is(linModsMPcell <- fitLMMs(objBG, features = getFeatures(objBG)[1:5],
        fixedVars = c("condition", "cellType"), pi = "midpoint"
    ), "list")
    # Account for cell as random effect
    expect_is(linModsEdgeCell <- fitLMMs(objBG, features = getFeatures(objBG)[1:5],
        fixedVars = "condition",
        randomVars = "root/cell", pi = "edge"
    ), "list")
    expect_is(linModsMidCellType <- fitLMMs(objBG, features = getFeatures(objBG)[1:5],
        fixedVars = c("condition", "cellType"),
        randomVars = "root/cell", pi = "midpoint"
    ), "list")
    expect_is(linModsNNCellType <- fitLMMs(objBG, features = getFeatures(objBG)[1:5],
        fixedVars = c("condition", "cellType"),
        pi = "nnCell", randomVars = "root"
    ), "list")
    expect_is(fitLMMsAll(objBG, fixedVars = c("condition", "cellType"),
        randomVars = "root", pis = c("nn", "nnCell"),
        features = getFeatures(objBG)[1:5]), "list")
    expect_is(resMat <- getResults(linModsMP, "Intercept"), "matrix")
    expect_is(resMatCond <- getResults(linModsEdge, "condition"), "matrix")
    expect_is(getResults(linModsMPcell, "cellType"), "matrix")
    expect_false(is.unsorted(getResults(linModsMP, "Intercept")[, "pVal"]))
})
test_that("Fitting linear mixed models throws errors where appropriate", {
    expect_error(fitLMMs(objBG,
        fixedVars = "treatment", randomVars = "fov",
        pi = "midpoint"
    ))
})
