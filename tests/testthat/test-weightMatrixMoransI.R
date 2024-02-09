context("Tests on Moran's I")
coordList = replicate(20,ppp("x" = runif(1), y = runif(1)), simplify = FALSE)
test_that("Weight matrix is correctly constructed", {
    wMat = buildMoransIWeightMat(coordList, numNNs = numNNs <- 8)
    expect_identical(unique(rowSums(wMat>0)), numNNs)
})
if(requireNamespace("ape")){
    m = 1e3;x = rnorm(m)
    coordMat = matrix(rnorm(2*m), ncol = 2)
    colnames(coordMat) = c("x", "y")
    nns = nnwhich(X = coordMat[, "x"], Y = coordMat[, "y"], k = seq.int(numNNs))
    W = sparseMatrix(i = rep(seq_len(nrow(coordMat)), numNNs), j = nns, x = 1/numNNs)
    #Prepare matrix as comes out of buildMoransIWeightMat
    test_that("Our implementation and that of the ape package match", {
        moranIape = ape::Moran.I(x, as.matrix(W), scaled = TRUE)
        expect_equal(moransI(x, W),
                         (moranIape$observed - moranIape$expected)/moranIape$sd
)
    })
}
