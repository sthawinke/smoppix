context("Test pim calculation spatrans package")
test_that("Calculating pims proceeds without errors", {
<<<<<<< HEAD
    expect_is(piEstsCSR <- estPims(hypFrame2,
        pis = c(
            "nn", "allDist",
            "nnPair", "allDistPair", "edge", "midpoint", "nnCell", "allDistCell",
            "nnPairCell", "allDistPairCell"
        ), null = "CSR",
        features = c("gene1", "gene2")
    ), "list")
    expect_is(piEstsBG <- estPims(hypFrame2,
        pis = c(
            "nn", "allDist", "nnPair", "allDistPair", "nnCell",
            "allDistCell", "nnPairCell", "allDistPairCell"
        ),
        null = "background", features = c("gene1", "gene2")
    ), "list")
=======
    expect_is(piEstsCSR <- estPims(hypFrame2, pis = c("nn", "allDist",
    "nnPair", "allDistPair", "edge", "midpoint", "nnCell", "allDistCell",
    "nnPairCell", "allDistPairCell"), null = "CSR",
    features = c("gene1", "gene2")), "list")
    expect_is(piEstsBG <- estPims(hypFrame2,
        pis = c("nn", "allDist", "nnPair", "allDistPair", "nnCell",
                "allDistCell", "nnPairCell", "allDistPairCell"),
        null = "background", features = c("gene1", "gene2")), "list")
>>>>>>> d2614d0d4780fc09e529e163278579595112568f
})
test_that("Calculating pims throws errors where appropriate", {
    expect_error(estPims(hypFrame2, pis = c("K-function")))
    expect_error(estPims(hypFrame, pis = c("edge")))
    expect_error(estPims(hypFrame2,
        pis = c("nn", "allDist", "nnPair", "allDistPair"),
        null = "background", features = c("gene200", "gene2")
    ))
})
