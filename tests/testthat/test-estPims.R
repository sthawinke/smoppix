context("Test pim calculation spatrans package")
test_that("Calculating pims proceeds without errors", {
    expect_s3_class(piEstsCSR <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint", "midpoint"),
                                    point = c(0.5, 0.5), null = "CSR"), "list")
    expect_s3_class(piEstsBG <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair"), null = "background"), "list")
    expect_message(piEstsBG2 <- estPims(hypFrame2, pis = c("edge", "fixedpoint", "midpoint"),
                                       point = c(0.5, 0.5), null = "background"))
})
test_that("Calculating pims throws errors where appropriate", {
    expect_error(estPims(hypFrame2, pis = c("K-function")))
    expect_error(estPims(hypFrame, pis = c("edge")))
})
