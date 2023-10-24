context("Test pim calculation spatrans package")
test_that("Calculating pims proceeds without errors", {
  expect_silent(piEsts <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint"),
                                    point = c(0.5, 0.5)))
})
test_that("Calculating pims throws errors where appropriate", {
    expect_error(estPims(hypFrame2, pis = c("K-function")))
    expect_error(estPims(hypFrame, pis = c("edge")))
})
