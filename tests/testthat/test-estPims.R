context("Test pim calculation spatrans package")

hypFrame <- buildHyperFrame(df, coordVars = c("x", "y"), designVar = "fov")
test_that("Calculating pims proceeds without errors", {
  expect_silent(piNN <- lapply(hypFrame$ppp, estPims, pis = c("nn", "allDist", "nnPair", "allDistPair")))
})
