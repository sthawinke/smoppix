context("Test PI calculation smoppix package")
test_that("Calculating pims proceeds without errors", {
    expect_is(piEstsCSR, "list")
    expect_is(piEstsBG, "list")
})
test_that("Calculating pims throws errors where appropriate", {
    expect_error(estPis(hypFrame2, pis = c("K-function")))
    expect_error(estPis(hypFrame, pis = c("edge")))
    expect_error(estPis(hypFrame2,
        pis = c("nn", "nnPair"), null = "background",
        features = c("gene200", "gene2")
    ))
    hypFrameUnsorted <- hypFrame2 # Unsort the x-coordinates
    coords(hypFrameUnsorted$ppp[[1]]) <- coords(hypFrameUnsorted$ppp[[1]])[sample(npoints(hypFrameUnsorted$ppp[[1]])), ]
    expect_error(estPis(hypFrameUnsorted, pis = "nn"))
})
