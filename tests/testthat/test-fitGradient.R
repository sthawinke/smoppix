context("Tests related to gradient fitting")
numPPPs = 10
hyp1 = hyperframe("ppp" = lapply(integer(numPPPs), function(i) runifpoint(1e2)))
Win = owin(poly = list(x = c(0, 0.5, 1, 0.5), y = c(0.5, 0, 0.5, 1)))
hyp2 = hyperframe("ppp" = lapply(integer(10), function(i) runifpoint(1e2, win = Win)))
test_that("fitGradient has correct return", {
    fg = fitGradient(hyp1, silent = TRUE, fixedForm = formula("ppp ~ id:x + id:y + id"),
                     fixedFormSimple = formula("ppp ~ id"), randomForm = NULL)
    expect_true(fg$pVal < 1 && fg$pVal > 0)
    expect_length(fg$coef, numPPPs*3)
    expect_s3_class(fitGradient(hyp1 ,silent = TRUE, fixedForm = formula("ppp ~ id:x + id:y + id"),
                                fixedFormReduced = formula("ppp ~ id"), randomForm = NULL,
                                returnModel = TRUE), "mppm")
})
test_that("fitGradient fails where appropriate", {
    expect_identical(fitGradient(runifpoint(5e1), silent = TRUE, fixedForm = formula("ppp ~ id:x + id:y + id"),
                                 fixedFormReduced = formula("ppp ~ id"), randomForm = NULL)$pVal, 1)
})
test_that("estGradients has correct return", {
    yangGrads <- estGradients(hypYang[seq_len(20),], features = feat <- getFeatures(hypYang)[seq_len(2)])
    expect_is(yangGrads, "list")
    expect_identical(names(yangGrads), feat)
})
test_that("estGradients throws errors where appropriate", {
    expect_error(estGradients(hypYang, gradients = "x"))
    expect_error(estGradients(hypEng[seq_len(10),], features = feat <- getFeatures(hypEng)[seq_len(2)],
                             fixedEffects = "treatment"))
})
test_that("estGradients works for cells as well", {
    engGrads <- estGradients(hypEng[seq_len(4),], features = feat <- getFeatures(hypEng)[seq_len(2)])
    expect_equal(names(engGrads), feat)
    expect_equal(names(engGrads[[1]]), c("overall", "cell"))
    expect_true(engGrads[[1]]$overall$pVal >= 0 & engGrads[[1]]$overall$pVal <= 1)
    engGradsCell <- estGradients(hypEng[seq_len(4),], features = feat <- getFeatures(hypEng)[seq_len(2)],
                             fixedEffects = "experiment")
    expect_true(engGradsCell[[1]]$cell$pVal >= 0 & engGradsCell[[1]]$cell$pVal <= 1)
    expect_length(engGradsCell[[1]]$cell$coef,
                  1+3*sum(vapply(hypEng[seq_len(4),]$ppp, function(x) length(unique(marks(x)$cell)), FUN.VALUE = 1)))
    #One extra parameter for experiment baseline
})

