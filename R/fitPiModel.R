fitPiModel = function(Formula, df, contrasts, Control, fixedPart){
    if (MM) {
        mod <- try(lmerTest::lmer(Formula, data = df, na.action = na.omit,
                                  weights = if (noWeight) {NULL} else {weight
                                  }, contrasts = contrasts, control = Control),
                   silent = TRUE)
    }
    # If no random variables or fit failed, go for linear model
    if (!MM || is(mod, "try-error")) {
        mod <- try(lm(formula(fixedPart), data = df, na.action = na.omit,
                      weights = if (noWeight) {
                          NULL
                      } else {
                          weight
                      }, contrasts = contrasts), silent = TRUE)
    }
    # If still fails, drop fixed effects
    if (is(mod, "try-error")) {
        mod <- lm(formula("pi - 0.5 ~ 1"), data = df, na.action = na.omit,
                  weights = if (noWeight) {
                      NULL
                  } else {
                      weight
                  })
    }
    return(mod)
}
