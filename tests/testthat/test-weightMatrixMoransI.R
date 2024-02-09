context("Tests on Moran's I")
coordList = replicate(20,list("x" = rnorm(1), y = rnorm(1)))
test_that("Weight matrix is correctly constructed", {
    wMat = buildWeightMat(coordList, numNNs = numNNs <- 8)
    expect_identical(unique(rowSums(wMat>0)), numNNs)
})
if(requireNamespace("ape")){
    m = 1e1;x = rnorm(m)
    W = matrix(runif(m^2), m , m)
    diag(W) = 0
    test_that("Our implementation and that of the ape package match", {
        moranIape = ape::Moran.I(x, W, scaled = TRUE)
        expect_identical(moranI(x, W),
                         (moranIape$observed - moranIape$expected)/moranIape$sd
)
    })
}
