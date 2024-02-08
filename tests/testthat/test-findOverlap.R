context("Overlap between windows")
library(spatstat.geom)
owins <- lapply(seq(0, 0.9, 0.1), function(i){
    owin(xrange = c(i, i + 0.25), yrange = c(0,1))
})
test_that("findOverlap really detects overlap between windows", {
    expect_error(findOverlap(owins))
    expect_is(idOverlap <- findOverlap(owins, returnIds = TRUE), "matrix")
})

