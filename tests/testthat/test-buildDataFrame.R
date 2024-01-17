context("Test building of dataframes of probabilistic indices for mixed model building")
test_that("Building data frames for mixed modelling proceeds without errors", {
    expect_silent(dfBG0 <- buildDataFrame(objBG, gene = "gene1", pi = "nn"))
    expect_s3_class(dfBG1 <- buildDataFrame(objBG, gene = "gene1", pi = "nn"), "data.frame")
    expect_s3_class(dfBG2 <- buildDataFrame(objBG, gene = "gene1--gene2", pi = "nnPair"), "data.frame")
    expect_s3_class(dfBG2b <- buildDataFrame(objBG, gene = c("gene1", "gene2"), pi = "nnPair"), "data.frame")
    expect_identical(dfBG2, dfBG2b)
    expect_silent(dfBG3 <- buildDataFrame(objBG, gene = "gene1", pi = "midpoint"))
    expect_silent(dfBG6 <- buildDataFrame(objBG, gene = "gene1", pi = "edge"))
    totPointsGene1 <- sum(vapply(objBG$hypFrame$ppp,
        FUN.VALUE = double(1),
        function(x) sum(marks(x, drop = TRUE)$gene == "gene1")
    ))
    # For fixed points, number of PI must match number of events
    expect_equal(totPointsGene1, nrow(dfBG6))
    expect_equal(sum(!is.na(dfBG0$pi)), sum(vapply(objBG$hypFrame$ppp,
        FUN.VALUE = logical(1),
        function(x) sum(marks(x, drop = TRUE)$gene == "gene1") > 1
    )))
    expect_s3_class(dfBG7 <- buildDataFrame(objBG, gene = "gene1", pi = "nnCell"), "data.frame")
    expect_s3_class(dfBG9 <- buildDataFrame(objBG, gene = "gene1--gene2", pi = "nnPairCell"), "data.frame")
    expect_silent(dfCSR1 <- buildDataFrame(objCSR, gene = "gene2", pi = "edge"))
    expect_silent(dfCSR3 <- buildDataFrame(objCSR, gene = "gene2", pi = "midpoint"))
    expect_s3_class(dfCSR5 <- buildDataFrame(objCSR, gene = "gene3", pi = "nn"), "data.frame")
    expect_s3_class(dfCSR6 <- buildDataFrame(objCSR, gene = "gene2--gene3", pi = "nnPair"), "data.frame")
    expect_s3_class(dfCSR8 <- buildDataFrame(objCSR, gene = "gene1", pi = "nnCell"), "data.frame")
    expect_s3_class(dfCSR9 <- buildDataFrame(objCSR, gene = "gene1--gene2", pi = "nnPairCell"), "data.frame")
    expect_true(all(is.na(dfCSR6$pi) | (dfCSR6$pi >= 0 & dfCSR6$pi <= 1)))
    expect_true(all(is.na(dfBG3$pi) | (dfBG3$pi >= 0 & dfBG3$pi <= 1)))
})
objCSR2 <- estPims(hypFrame2,
    pis = c("nn"), features = c("gene1", "gene2"),
    null = "CSR"
)
test_that("Building data frames throws errors where appropriate", {
    expect_error(buildDataFrame(objCSR, gene = c("gene1", "gene2"), pi = "nn"))
    expect_error(buildDataFrame(objCSR, gene = "gene1", pi = "Kest"))
    expect_error(buildDataFrame(objBG, gene = "protein1", pi = "nn"))
    expect_error(buildDataFrame(objCSR2, gene = "gene1--gene2", pi = "nnPair"))
    expect_error(buildDataFrame(objBG, gene = c("gene1", "gene2"), i = "nn"))
    expect_error(buildDataFrame(objBG, gene = "gene1_gene2", i = "nn"))
    expect_error(buildDataFrame(objBG, gene = "gene1", i = "nnPair"))
    expect_error(buildDataFrame(objBG, gene = c("gene1", "gene2", "gene3"), pi = "nnPair", hypFrame = hypFrame2))
})
