context("Plotting the weight function")
data(Yang)
hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
designVar = c("day", "root", "section"))
yangPims = estPims(hypYang, pis = c("nn", "nnPair"),
                   features = attr(hypYang, "features")[1:20])
#First Build the weight function
wf <- buildWeightFunction(yangPims, pi = "nn", hypFrame = hypYang,
                          designVars = c("day", "root"))
wfPair <- buildWeightFunction(yangPims, pi = "nnPair", hypFrame = hypYang,
designVars = c("day", "root"))
test_that("Plotting the weight function works", {
    expect_silent(plotWf(wf))
    expect_is(tmpPlot <- plotWf(wfPair), "ggplot")
})
test_that("Plotting the weight function application throws errors where appropriate", {
    expect_error(plotWf(piEstsBG))
})
