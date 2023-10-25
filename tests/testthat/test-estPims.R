context("Test pim calculation spatrans package")
test_that("Calculating pims proceeds without errors", {
    expect_silent(piEstsCSR <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint", "midpoint"),
                                    point = c(0.5, 0.5), null = "CSR"))
    expect_silent(piEstsBG <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair"), null = "background"))
    expect_message(piEstsCSR2 <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint", "midpoint"),
                                       point = c(0.5, 0.5), null = "background"))
})
test_that("Calculating pims throws errors where appropriate", {
    expect_error(estPims(hypFrame2, pis = c("K-function")))
    expect_error(estPims(hypFrame, pis = c("edge")))
})
