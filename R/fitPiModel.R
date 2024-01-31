#' @importFrom lmerTest lmer
#' @importFrom stats na.omit lm
fitPiModel = function(Formula, dff, contrasts, Control, MM){
    if (MM) {
        mod <- try(lmerTest::lmer(Formula, data = dff, na.action = na.omit,
                                  weights = weight, contrasts = contrasts, control = Control),
                   silent = TRUE)
    }
    # If no random variables or fit failed, go for linear model
    if (!MM || is(mod, "try-error")) {
        mod <- try(lm(getFixedPart(Formula), data = dff, na.action = na.omit,
                      weights = weight, contrasts = contrasts), silent = TRUE)
    }
    return(mod)
}
