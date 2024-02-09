context("Fitting models works in any setting")
n <- 100
library(lmerTest)
library(lme4)
df <- data.frame(outcome = rnorm(n), fixed = sample(c(TRUE, FALSE), n, replace = TRUE), random = sample(paste0(
    "subject",
    seq_len(3)
), n, replace = TRUE))
weight <- runif(n)
weight <- weight / sum(weight)
Control <- lmerControl()
contrasts <- list(fixed = "contr.sum")
test_that("Mixed model fitting works", {
    expect_s4_class(fitPiModel("outcome~fixed + (1|random)", df, contrasts, Control, MM = TRUE), "lmerModLmerTest")
    expect_s4_class(fitPiModel("outcome~fixed + (1|random)", df, contrasts, Control, MM = TRUE, Weight = weight), "lmerModLmerTest")
})
test_that("Fixed model fitting works", {
    expect_s3_class(fitPiModel("outcome~fixed", df, contrasts, Control, MM = FALSE), "lm")
    expect_s3_class(fitPiModel("outcome~fixed", df, contrasts, Control, MM = FALSE, Weight = weight), "lm")
})
