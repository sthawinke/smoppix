context("Plotting of top results")
linModsMP <- fitLMMs(objBG, returnModels = TRUE, features = getFeatures(objBG)[1:5],
    fixedVars = "condition", pi = "centroid"
)
linModsNNint <- fitLMMs(yangPims, fixedVars = "day", randomVars = "root", features = getFeatures(yangPims)[15:20])
test_that("Top results plotting proceeds without errors", {
    expect_silent(plotTopResults(results = linModsNNint, hypYang, "nn"))
    expect_silent(plotTopResults(
        results = linModsNNint, hypYang, "nnPair", effect = "day",
        sigLevel = 0.99, effectParameter = "day0"
    ))
    expect_silent(plotTopResults(
        results = linModsNNint, hypYang, "nnPair", effect = "day",
        sigLevel = 0.99, effectParameter = "day3"
    ))
    expect_silent(plotTopResults(
        results = linModsNNint, hypYang, "nnPair", what = "anti",
        sigLevel = 0.4
    ))
})
test_that("Fitting linear mixed models throws errors where appropriate", {
    expect_error(getResults(results = linModsNNint, hypYang, "centroid"))
    expect_error(getResults(results = linModsNNint, hypYang, "nn", effect = "treatment"))
    expect_error(plotTopResults(results = linModsNNint, hypYang, "centroid"))
    expect_error(plotTopResults(results = linModsNNint, hypYang, "nn", effect = "treatment"))
    expect_error(plotTopResults(results = linModsMP, hypFrame, "centroid", sigLevel = 1e-05))
    expect_error(plotTopResults(results = linModsNNint, hypYang, "nn", what = "repulsion"))
    expect_error(plotTopResults(results = linModsNNint, hypYang, "nn", what = "regular", effect = "day"))
})
