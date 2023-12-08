#' Fit linear mixed nodels for all features of a pimRes object
#'
#' @inheritParams buildDfMM
#' @param fixedVars,randomVars Names of fixed and random variables
#' @param verbose A boolean, should the formula be printed?
#' @param returnModels a boolean: should the full models be returned?
#' Otherzise only summary statistics are returned
#' @param Formula A formula; if not supplied it will be constructed
#' from the fixed and random variables
#'
#' @return A fitted linear mixed model of class 'lmerTest'
#' @export
#' @importFrom lmerTest lmer
#' @importFrom stats formula na.omit anova
#' @importFrom lme4 lmerControl .makeCC isSingular
#'
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'   coordVars = c("x", "y"),
#'   imageVars = c("day", "root", "section")
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPims(hypYang, pis = "nn", features = attr(hypYang, "features")[seq_len(20)])
#' # First build the weight function
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' lmmModels <- fitLMMs(yangObj, fixedVars = "day", randomVars = "root")
fitLMMs <- function(resList, pi, fixedVars = NULL, randomVars = NULL, verbose = TRUE,
                    returnModels = FALSE, Formula = NULL) {
  pi <- match.arg(pi, choices = c(
    "nn", "allDist", "nnPair", "allDistPair",
    "edge", "midpoint", "nnCell", "allDistCell",
    "nnPairCell", "allDistPairCell"
  ))
  designVars <- c(fixedVars, randomVars)
  stopifnot(all(designVars %in% resList$designVars))
  if (is.null(Formula)) {
    Formula <- formula(formChar <- paste(
      "pi - 0.5 ~ 1",
      paste(
        if (!is.null(fixedVars)) {
          paste("+", paste(fixedVars, sep = "+"))
        }, if (!is.null(randomVars)) {
          paste("+", paste0("(1|", randomVars, ")", collapse = "+"))
        }
      )
    ))
  }
  MM <- any(grepl("|", formChar))
  if (verbose) {
    message("Fitted formula:\n", formChar)
  }
  Control <- lmerControl( ## convergence checking options
    check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL),
    check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),
    check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6)
  )
  if (is.null(fixedVars)) {
    contrasts <- NULL
  } else {
    contrasts <- lapply(fixedVars, function(x) "contr.sum")
    names(contrasts) <- fixedVars
  }
  Features <- if (grepl("Pair", pi)) {
    makePairs(attr(resList$hypFrame, "featuresEst"))
  } else {
    attr(resList$hypFrame, "featuresEst")
  }
  models <- lapply(Features, function(gene) {
    df <- buildDfMM(resList, gene = gene, pi = pi)
    if (sum(!is.na(df$pi)) < 3) {
      return(NULL)
    }
    if (MM) {
      mod <- try(lmerTest::lmer(Formula,
        data = df, na.action = na.omit, weights = weight,
        contrasts = contrasts, control = Control
      ), silent = TRUE)
    }
    # If no random variables or fit failed, go for linear model
    if (!MM || is(mod, "try-error")) {
      mod <- lm(Formula,
        data = df, na.action = na.omit, weights = weight,
        contrasts = contrasts
      )
    }
    return(mod)
  })
  results <- extractResults(models, fixedVars)
  # Effect size, standard error, p-value and adjusted p-value per mixed effect
  if (!returnModels) {
    return(list("results" = results, "resList" = resList))
  } else {
    return(list("results" = results, "resList" = resList, "models" = models))
  }
}
