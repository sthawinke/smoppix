context("Test building of dataframes of probablistic indices for mixed model building")
piEstsCSR <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint", "midpoint"),
                     point = c(0.5, 0.5), null = "CSR")
piEstsBG <- estPims(hypFrame2, pis = c("edge", "fixedpoint", "midpoint"),
                      point = c(0.5, 0.5), null = "background")
test_that("Building data frames for mixed modelling proceeds without errors", {
    expect_s3_class(dfCSR1 <- buildDfMM(piEstsCSR, gene = "gene1", pi = "nn"), "data.frame")
    expect_s3_class(dfCSR2 <- buildDfMM(piEstsCSR, gene = "gene1", pi = "nnPair"), "data.frame")
    expect_s3_class(dfCSR3 <- buildDfMM(piEstsCSR, gene = "gene1", pi = "allDist"), "data.frame")
    expect_s3_class(dfCSR4 <- buildDfMM(piEstsCSR, gene = "gene1", pi = "allDistPair"), "data.frame")
    expect_s3_class(dfBG1 <- buildDfMM(piEstsBG, gene = "gene1", pi = "edge"), "data.frame")
    expect_s3_class(dfBG2 <- buildDfMM(piEstsBG, gene = "gene3", pi = "fixedpoint"), "data.frame")
    expect_s3_class(dfBG3 <- buildDfMM(piEstsBG, gene = "gene2", pi = "midpoint"), "data.frame")
})
test_that("Calculating pims throws errors where appropriate", {
    expect_error(buildDfMM(piEstsCSR, gene = c("gene1", "gene2"), pi = "nn"))
    expect_error(buildDfMM(piEstsCSR, gene = "gene1", pi = "Kest"))
    expect_error(buildDfMM(piEstsBG, gene = "gene1", pi = "nn"))
})
