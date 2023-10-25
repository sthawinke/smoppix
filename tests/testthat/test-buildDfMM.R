context("Test building of dataframes of probablistic indices for mixed model building")
piEstsCSR <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint", "midpoint"),
                     point = c(0.5, 0.5), null = "CSR")
piEstsBG <- estPims(hypFrame2, pis = c("edge", "fixedpoint", "midpoint"),
                      point = c(0.5, 0.5), null = "background")
test_that("Building data frames fro mixed modelling proceeds without errors", {
    expect_s3_class(dfCSR <- buildDfMM(piEstsCSR, gene = "gene1", pi = "nn"), "data.frame")
    expect_s3_class(dfBG <- buildDfMM(piEstsBG, gene = "gene1", pi = "nn"), "data.frame")
})
test_that("Calculating pims throws errors where appropriate", {
    expect_error(buildDfMM(piEstsCSR, gene = c("gene1", "gene2"), pi = "nn"))
    expect_error(buildDfMM(piEstsCSR, gene = "gene1", pi = "Kest"))
    expect_error(buildDfMM(piEstsBG, gene = "gene1", pi = "nn"))
})
