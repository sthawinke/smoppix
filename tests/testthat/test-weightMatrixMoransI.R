context("Tests on Moran's I")
library(Matrix)
coordMat <- matrix(rnorm(40), ncol = 2)
numNNs <- 8
test_that("Weight matrix is correctly constructed", {
    wMat <- buildMoransIWeightMat(coordMat, numNNs = numNNs)
    expect_equal(unique(Matrix::rowSums(wMat > 0)), numNNs)
})
if (requireNamespace("ape")) {
    m <- 1000
    x <- rnorm(m)
    coordMat <- matrix(rnorm(2 * m), ncol = 2)
    colnames(coordMat) <- c("x", "y")
    nns <- nnwhich(X = coordMat[, "x"], Y = coordMat[, "y"], k = seq.int(numNNs))
    W <- sparseMatrix(i = rep(seq_len(nrow(coordMat)), numNNs), j = nns, x = 1/numNNs)
    # Prepare matrix as comes out of buildMoransIWeightMat
    test_that("Our implementation and that of the ape package match", {
        moranIape <- ape::Moran.I(x, as.matrix(W), scaled = TRUE)
        expect_equal(unname(moransI(x, W)["MoransI"]), (moranIape$observed - moranIape$expected))
    })
}
