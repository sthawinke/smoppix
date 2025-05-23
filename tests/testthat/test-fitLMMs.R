context("Large scale mixed model fitting")
test_that("Fitting linear mixed models proceeds without errors", {
    expect_is(
        linModsNNint <- fitLMMs(yangPims, pi = "nn", features = getFeatures(yangPims)[15:18]),
        "list"
    )
    expect_is(
        fitLMMs(yangPims, pi = "nn", verbose = FALSE, features = getFeatures(yangPims)[15:18]),
        "list"
    )
    # Supply your own formula
    expect_is(fitLMMs(yangPims,
        features = getFeatures(yangPims)[15:18], pi = "nn",
        Formula = "pi - 0.5 ~ day +1|root"
    ), "list")
    expect_is(linModsNNPairint <- fitLMMs(yangPims,
        features = getFeatures(yangPims)[15:18], pi = "nnPair"
    ), "list")
    expect_is(linModsNN <- fitLMMs(yangPims,
        features = getFeatures(yangPims)[15:18],
        fixedVars = "day", pi = "nn"
    ), "list")
    expect_is(linModsNNPair <- fitLMMs(yangPims, 
                                       features = getFeatures(yangPims)[15:18],
        fixedVars = "day", pi = "nnPair"
    ), "list")
    expect_is(linMModsNN <- fitLMMs(yangPims,
                                    features = getFeatures(yangPims)[15:18],
        fixedVars = "day", randomVars = "root", pi = "nn"
    ), "list")
    expect_is(linMModsNNPair <- fitLMMs(yangPims,
                                        features = getFeatures(yangPims)[15:18],
        fixedVars = "day", randomVars = "root", pi = "nnPair"
    ), "list")
    # Returning the models
    expect_is(linModsNNfull <- fitLMMs(yangPims,features = getFeatures(yangPims)[15:18],
        fixedVars = "day", pi = "nn", returnModels = TRUE
    ), "list")
    expect_is(
        linMModsNNfull <- fitLMMs(yangPims,
            features = getFeatures(yangPims)[15:18],
            fixedVars = "day", randomVars = "root", pi = "nn", returnModels = TRUE
        ),
        "list"
    )
    expect_s3_class(linModsNNfull[["nn"]]$models[[getFeatures(yangPims)[[15]]]]$piMod, "lm")
    expect_s4_class(linMModsNNfull[["nn"]]$models[[getFeatures(yangPims)[[15]]]]$piMod, "lmerModLmerTest")
    expect_is(linModsMP <- fitLMMs(objBG,
        returnModels = TRUE, features = getFeatures(objBG)[1:5],
        fixedVars = "condition", pi = "centroid", addMoransI = TRUE
    ), "list")
    expect_is(getResults(linModsMP, "centroid", "Intercept", moransI = TRUE), "matrix")
    expect_s3_class(linModsMP[["centroid"]]$models[[1]]$moranMod, "lm")
    expect_s4_class(linModsMP[["centroid"]]$models[[1]]$piMod, "lmerModLmerTest")
    expect_is(linModsEdge <- fitLMMs(objBG,
        features = getFeatures(objBG)[1:5],
        fixedVars = c("condition", "age"), pi = "edge"
    ), "list")
    expect_identical(ncol(getResults(linModsEdge, "edge", "age")), 3L)
    # Including cell (type) either as fixed or random effect E.g. test for
    # differences between cell types
    expect_is(linModsMPcell <- fitLMMs(objBG,
        features = getFeatures(objBG)[1:5],
        fixedVars = c("condition", "cellType"), pi = "centroid"
    ), "list")
    # Account for cell as random effect
    expect_is(linModsEdgeCell <- fitLMMs(objBG,
        features = getFeatures(objBG)[1:5],
        fixedVars = "condition", randomVars = "image/cell", pi = "edge"
    ), "list")
    expect_is(linModsMidCellType <- fitLMMs(objBG,
        features = getFeatures(objBG)[1:5],
        addMoransI = TRUE, fixedVars = c("condition", "cellType"), pi = "centroid",
        returnModels = TRUE,
    ), "list")
    expect_s4_class(linModsMidCellType[["centroid"]]$models[[1]]$piMod, "lmerModLmerTest")
    expect_s3_class(linModsMidCellType[["centroid"]]$models[[1]]$moranMod, "lm")
    expect_is(linModsNNCellType <- fitLMMs(objBG,
        features = getFeatures(objBG)[1:3],
        fixedVars = c("condition", "cellType"), pi = "nnCell",
    ), "list")
    expect_warning(fitLMMs(objBG, fixedVars = c("condition", "cellType"), pis = c(
        "nn",
        "nnCell"
    ), features = getFeatures(objBG)[1:3]))
    expect_is(resMat <- getResults(linModsMP, "centroid", "Intercept"), "matrix")
    expect_is(resMatCond <- getResults(linModsEdge, "edge", "condition"), "matrix")
    expect_is(getResults(linModsMPcell, "centroid", "cellType"), "matrix")
    expect_false(is.unsorted(getResults(linModsMP, "centroid", "Intercept")[, "pVal"]))
})
test_that("Fitting linear mixed models throws errors where appropriate", {
    expect_error(fitLMMs(objBG, fixedVars = "treatment", randomVars = "fov", pi = "centroid"))
})
