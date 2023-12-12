context("Large scale mixed model fitting")
test_that("Fitting linear mixed models proceeds without errors", {
  expect_message(linModsNNint <- fitLMMs(yangObj, pi = "nn"))
  expect_silent(fitLMMs(yangObj, pi = "nn", verbose = FALSE))
  # Supply your own formula
  expect_message(fitLMMs(yangObj,
    pi = "nn",
    Formula = "pi - 0.5 ~ day +1|root"
  ))
  expect_message(linModsNNPairint <- fitLMMs(yangObj, pi = "nnPair"))
  expect_message(linModsNN <- fitLMMs(yangObj, fixedVars = "day", pi = "nn"))
  expect_message(linModsNNPair <- fitLMMs(yangObj,
    fixedVars = "day",
    pi = "nnPair"
  ))
  expect_message(linMModsNN <- fitLMMs(yangObj,
    fixedVars = "day",
    randomVars = "root", pi = "nn"
  ))
  expect_message(linMModsNNPair <- fitLMMs(yangObj,
    fixedVars = "day",
    randomVars = "root", pi = "nnPair"
  ))
  # Returning the models
  expect_message(linModsNNfull <- fitLMMs(yangObj,
    fixedVars = "day",
    pi = "nn", returnModels = TRUE
  ))
  expect_message(linMModsNNfull <- fitLMMs(yangObj,
    fixedVars = "day",
    randomVars = "root", pi = "nn", returnModels = TRUE
  ))
  expect_s3_class(linModsNNfull$models[[1]], "lm")
  expect_s4_class(linMModsNNfull$models[[1]], "lmerModLmerTest")
  expect_message(linModsMP <- fitLMMs(objBG,
    fixedVars = "condition",
    pi = "midpoint"
  ))
  expect_message(linModsEdge <- fitLMMs(objBG,
    fixedVars = "condition",
    pi = "edge"
  ))
  # Including cell (type) either as fixed or random effect
  # E.g. test for differences between cell types
  expect_message(linModsMPcell <- fitLMMs(objBG,
    fixedVars = c("condition", "cellType"), pi = "midpoint"
  ))
  # Check nesting for lmerTest! FIX ME!
  # Account for cell as random effect
  expect_message(linModsEdgeCell <- fitLMMs(objBG,
    fixedVars = "condition",
    randomVars = "cell", pi = "edge"
  ))
  expect_message(linModsMidCellType <- fitLMMs(objBG,
    fixedVars = c("condition", "cellType"),
    randomVars = "cell", pi = "midpoint"
  ))
  expect_message(linModsNNCellType <- fitLMMs(objBG,
    fixedVars = c("condition", "cellType"),
    pi = "nnCell"
  ))
  expect_is(resMat <- getResults(linModsMP, "Intercept"), "matrix")
  expect_is(resMatCond <- getResults(linModsEdge, "condition"), "matrix")
  expect_is(getResults(linModsMPcell, "cellType"), "matrix")
  expect_false(is.unsorted(getResults(linModsMP, "Intercept")[, "pVal"]))
})
test_that("Fitting linear mixed models throws errors where appropriate", {
  expect_error(fitLMMs(objBG,
    fixedVars = "condition", randomVars = "fov",
    pi = "midpoint"
  ))
})
