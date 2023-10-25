context("Test building of dataframes of probablistic indices for mixed model building")
piEstsCSR <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint", "midpoint"),
                     point = c(0.5, 0.5), null = "CSR")
piEstsCSR2 <- estPims(hypFrame2, pis = c("edge", "fixedpoint", "midpoint"),
                      point = c(0.5, 0.5), null = "background")
test_that("Building data frames fro mixed modelling proceeds without errors", {
    expect_silent(dfCSR <- buildDfMM(piEstsCSR, gene = "gene1", pi = "nn"))
})
test_that("Calculating pims throws errors where appropriate", {
    expect_error(estPims(hypFrame2, pis = c("K-function")))
    expect_error(estPims(hypFrame, pis = c("edge")))
})
