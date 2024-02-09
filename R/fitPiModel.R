#' Fit a linear model for an individual gene and PI combination
#'
#' @param dff The dataframe
#' @param contrasts The contrasts to be used, see \link{model.matrix}
#' @param Control Control parameters
#' @param MM A boolean, should a mixed model be tried
#' @param Weight A weight variable
#' @inheritParams fitLMMsSingle
#' @return A fitted model
#'
#' @importFrom lmerTest lmer
#' @importFrom stats na.omit lm as.formula
#' @importFrom lme4 nobars
#' @seealso \link{fitLMMsSingle}
fitPiModel <- function(Formula, dff, contrasts, Control, MM, Weight = NULL) {
    ff <- as.formula(Formula)
    environment(ff) <- environment() # Make sure formula sees the weights, see
    # https://stackoverflow.com/questions/61164404/call-to-weight-in-lm-within-function-doesnt-evaluate-properly
    if (MM) {
        mod <- try(lmerTest::lmer(ff, data = dff, na.action = na.omit, weights = Weight, contrasts = contrasts, control = Control),
            silent = TRUE
        )
    }
    # If no random variables or fit failed, go for linear model
    if (!MM || is(mod, "try-error")) {
        mod <- try(lm(nobars(ff), data = dff, na.action = na.omit, weights = Weight, contrasts = contrasts), silent = TRUE)
    }
    return(mod)
}
