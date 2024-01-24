context("Plotting of top results")
test_that("Top results plotting proceeds without errors", {
    linModsNNint <- fitLMMs(yangPims, features = getFeatures(yangPims)[seq_len(10)])
    expect_silent(plotTopResults(results = linModsNNint, hypYang, "nn"))
    expect_silent(plotTopResults(results = linModsNNint, hypYang, "nnPair"))
    expect_silent(plotTopResults(results = linModsNNint, hypYang, "nnPair", smallPI = FALSE))
})
test_that("Fitting linear mixed models throws errors where appropriate", {
    expect_error(plotTopResults(linModsNNint, hypYang, "centroid"))
})
