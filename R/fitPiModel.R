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
        mod <- try(lmerTest::lmer(ff, data = dff, na.action = na.omit, 
                                  weights = Weight, contrasts = contrasts, 
                                  control = Control), silent = TRUE)
    }
    # If no random variables or fit failed, go for linear model
    if (!MM || is(mod, "try-error")) {
        mod <- try(lm(nobars(ff),
            data = dff, na.action = na.omit, weights = Weight,
            contrasts = contrasts
        ), silent = TRUE)
    }
    return(mod)
}
#' Take an existing frame, add outcome and weight and fit lmer model
#'
#' @param ff The prepared frame
#' @param y outcome vector
#' @param weights weights vector
#'
#' @returns A fitted lmer model
#' @importFrom lme4 mkLmerDevfun optimizeLmer mkMerMod
fitSingleLmmModel <- function(ff, y, Control, weights = NULL) {
  fr <- ff$fr                    # this is a data.frame (model frame)
  ## Use model-frame column names used by stats::model.frame
  weights <- if (is.null(weights)) rep(1, length(y)) else weights
  #Set weights to zero for NA y's. This is more efficient than dropping rows, as we can keep using the same design matrix
  weights[naId <- is.na(y)] = 0
  y[naId] = 0 #Set to arbitrary number to avoid errors, will be downweighted anyway
  fr$`(weights)` = weights
  fr[["pi - 0.5"]] <- y                 # replace response
  
  devfun <- mkLmerDevfun(fr, ff$X, ff$reTrms, control = Control)
  opt    <- optimizeLmer(devfun, control = Control)
  mkMerMod(environment(devfun), opt, ff$reTrms, fr)
  # Switch to fixed effects model when failed = To DO!
}