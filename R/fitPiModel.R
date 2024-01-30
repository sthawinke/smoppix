#' @importFrom lmerTest lmer
#' @importFrom stats na.omit lm
fitPiModel = function(Formula, df, contrasts, Control, MM){
    if (MM) {
        mod <- try(lmerTest::lmer(Formula, data = df, na.action = na.omit,
                                  weights = if (is.null(df$weight)) {NULL} else {weight
                                  }, contrasts = contrasts, control = Control),
                   silent = TRUE)
    }
    # If no random variables or fit failed, go for linear model
    if (!MM || is(mod, "try-error")) {
        mod <- try(lm(getFixedPart(Formula), data = df, na.action = na.omit,
                      weights = if (is.null(df$weight)) {
                          NULL
                      } else {
                          weight
                      }, contrasts = contrasts), silent = TRUE)
    }
    return(mod)
}
