context("Test pim calculation spatrans package")
test_that("Calculating pims proceeds without errors", {
    expect_s3_class(piEstsCSR <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint", "midpoint"),
                                    point = c(0.5, 0.5), null = "CSR", features = c("gene1", "gene2")), "list")
    expect_s3_class(piEstsBG <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair"),
                                        null = "background", features = c("gene1", "gene2")), "list")
})
test_that("Calculating pims throws errors where appropriate", {
    expect_error(estPims(hypFrame2, pis = c("K-function")))
    expect_error(estPims(hypFrame, pis = c("edge")))
    expect_error(estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair"),
                                        null = "background", features = c("gene200", "gene2")), "list")
})
