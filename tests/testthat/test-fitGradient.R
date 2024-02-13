context("Tests related to gradient fitting")
ppp1 = runifpoint(1e2)
ppp2 = ppp1
ppp2$window = owin(poly = list(x = c(0, 0.5, 1, 0.5), y = c(0.5, 0, 0.5, 1)))
test_that("fitGradient has correct return", {
    expect_is(fitGradient(ppp1), "double")
    expect_s3_class(fitGradient(ppp1, returnModel = TRUE), "ppm")
})
test_that("fitGradient throws errors where appropriate", {
    expect_error(fitGradient(coords(ppp1)))
})
