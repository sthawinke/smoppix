context("Adding and lotting the weight function")
data(Yang)
<<<<<<< HEAD
hypYang <- buildHyperFrame(Yang,
    coordVars = c("x", "y"),
    imageVars = c("day", "root", "section")
)
yangPims <- estPims(hypYang,
    pis = c("nn", "nnPair"),
    features = attr(hypYang, "features")[1:20]
)
# First build the weight function
=======
hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
imageVars = c("day", "root", "section"))
yangPims = estPims(hypYang, pis = c("nn", "nnPair"),
                   features = attr(hypYang, "features")[1:20])
#First build the weight function
>>>>>>> d2614d0d4780fc09e529e163278579595112568f
yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
test_that("Plotting the weight function works", {
    expect_silent(plotWf(yangObj, "nn"))
    expect_is(tmpPlot <- plotWf(yangObj, "nnPair"), "ggplot")
})
test_that("Plotting the weight function application throws errors where appropriate", {
    expect_error(plotWf(piEstsBG))
})
