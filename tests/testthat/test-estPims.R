context("Test pim calculation spatrans package")
test_that("Calculating pims proceeds without errors", {
    expect_is(piEstsCSR, "list")
    expect_is(piEstsBG, "list")
})
test_that("Calculating pims throws errors where appropriate", {
    expect_error(estPims(hypFrame2, pis = c("K-function")))
    expect_error(estPims(hypFrame, pis = c("edge")))
    expect_error(estPims(hypFrame2,
        pis = c("nn", "nnPair"),
        null = "background", features = c("gene200", "gene2")
    ))
})
