context("Adding and lotting the weight function")
data(Yang)
hypYang <- buildHyperFrame(Yang,
  coordVars = c("x", "y"),
  imageVars = c("day", "root", "section")
)
yangPims <- estPims(hypYang,
  pis = c("nn", "nnPair"),
  features = attr(hypYang, "features")[1:20]
)
# First build the weight function
yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
test_that("Plotting the weight function works", {
  expect_silent(plotWf(yangObj, "nn"))
  expect_is(tmpPlot <- plotWf(yangObj, "nnPair"), "ggplot")
})
test_that("Plotting the weight function application throws errors where appropriate", {
  expect_error(plotWf(piEstsBG))
})
