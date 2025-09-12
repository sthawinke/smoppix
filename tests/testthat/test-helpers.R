context("Unit tests for small auxiliary function")
test_that("sund yields correct result", {
    expect_identical(sund("gene1--gene2"), c("gene1", "gene2"))
})
test_that("Make gene pairs works as expected", {
    expect_identical(makePairs(paste0("gene", seq_len(3))), c(
        "gene1--gene2", "gene1--gene3",
        "gene2--gene3"
    ))
})
test_that("getGp does correct extraction", {
    geneMat <- matrix(rnorm(6), 3, 3)
    colnames(geneMat) <- makePairs(paste0("gene", seq_len(3)))
    expect_identical(getGp(geneMat, "gene1--gene2"), getGp(geneMat, "gene2--gene1"))
    expect_is(getGp(geneMat, "gene1--gene2", drop = FALSE), "matrix")
})
library(lme4)
test_that("Design variables are well constructed", {
    expect_identical(makeDesignVar(data.frame(a = 1:3, b = LETTERS[1:3]), c(
        "a",
        "b"
    )), c("1_A", "2_B", "3_C"))
})
test_that("Coordinate extraction works", {
    expect_is(mat <- getCoordsMat(hypYang$ppp[[1]]), "matrix")
    expect_is(mat <- getCoordsMat(mat), "matrix")
    expect_is(mat <- getCoordsMat(as.data.frame(mat)), "matrix")
})
test_that("Feature pairs are correctly sorted", {
  genePairs <- c("geneB--geneC", "geneC--geneA", "gene2--gene1")
  genePairsSorted <- c("geneB--geneC", "geneA--geneC", "gene1--gene2")
  names(genePairsSorted) = genePairsSorted
  expect_identical(sortGp(genePairs), genePairsSorted)
})
