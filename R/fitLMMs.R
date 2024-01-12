#' Fit linear mixed nodels for all features of a pimRes object
#'
#' @inheritParams buildDataFrame
#' @param fixedVars Names of fixed effects
#' @param randomVars Names of random variables, possibly in a single vector to
#' reflect nesting structure, see details.
#' @param verbose A boolean, should the formula be printed?
#' @param returnModels a boolean: should the full models be returned?
#' Otherwise only summary statistics are returned
#' @param Formula A formula; if not supplied it will be constructed
#' from the fixed and random variables
#' @param randomNested A boolean, indicating if random effects are nested within
#' point patterns. See details.
#' @param features The features for which to fit linear mixed models.
#' Defaults to all features in the object
#' @details Genes or gene pairs with insufficient observations will be silently
#' omitted. When randomVars is provided as a vector, independent random
#' intercepts are fitted for them by default. Providing them separated by '\' or
#' ':' as in the lmer formulas is also allowed to reflect nesting structure.
#'
#' It is by default assumed that random effects are nested within the point
#'  patterns. This means for instance that cells with the same name but from
#'  different point patterns are assigned to different random effects. Set
#'  'randomNested' to FALSE to override this behaviour.
#'
#' @return A list of test results, if requested also the linear models are returned
#' @importFrom lmerTest lmer
#' @importFrom stats formula na.omit anova lm
#' @importFrom lme4 lmerControl .makeCC isSingular
#' @importFrom methods is
#' @seealso \link{buildDataFrame},\link{getResults}
fitLMMsSingle <- function(obj, pi, fixedVars = NULL, randomVars = NULL,
                          verbose = TRUE, returnModels = FALSE, Formula = NULL,
                          randomNested = TRUE, features = getFeatures(obj)) {
    pi <- match.arg(pi, choices = c(
        "nn", "nnPair", "edge", "midpoint",
        "nnCell", "nnPairCell"
    ))
    noWeight <- pi %in% c("edge", "midpoint")
    if (noWeight && !("cell" %in% randomVars)) {
        randomVars <- c(randomVars, "cell")
    }
    # For independent distances, no weights are needed
    randomVarsSplit <- if (!is.null(randomVars)) {
        unlist(lapply(c("/", ":"), function(Split) {
            tmp <- strsplit(randomVars, Split)
            if (length(tmp) > 2) {
                tmp
            }
        }))
    }
    designVars <- c(fixedVars, randomVarsSplit)
    if (any(id <- !(designVars %in% getDesignVars(obj)))) {
        stop("Design variables", designVars[id], "not found in object.")
    }
    # Allow cell as design variable, both fixed and random
    if (is.null(Formula)) {
        fixedPart <- paste("pi - 0.5 ~ 1", paste(if (!is.null(fixedVars)) {
            paste("+", paste(fixedVars, collapse = " + "))
        }))
        Formula <- formula(formChar <- paste(fixedPart, if (!is.null(randomVars)) {
            paste("+", paste0("(1|", randomVars, ")", collapse = " + "))
        }))
    } else {
        formChar <- paste(as.character(Formula), collapse = "")
    }
    MM <- any(grepl("\\|", formChar))
    if (verbose) {
        message("Fitted formula for pi ", pi, ":\n", formChar)
    }
    Control <- lmerControl(check.conv.grad = .makeCC("ignore", tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(
        action = "ignore",
        tol = formals(isSingular)$tol
    ), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06))
    # convergence checking options
    if (is.null(fixedVars)) {
        contrasts <- NULL
    } else {
        names(fixedVars) <- fixedVars
        contrasts <- lapply(fixedVars, function(x) named.contr.sum)
    }
    tmp <- if (grepl("Pair", pi)) {
        features <- makePairs(features)
        vapply(features, FUN.VALUE = double(nrow(obj$hypFrame)), function(gene) {
            vapply(obj$hypFrame$tabObs, FUN.VALUE = double(1), function(x) all(x[sund(gene)] >= 1))
        })
    } else {
        # First check if gene is present at all
        vapply(features, FUN.VALUE = double(nrow(obj$hypFrame)), function(gene) {
            vapply(obj$hypFrame$tabObs, FUN.VALUE = double(1), function(x) x[gene])
        })
    }
    featIds <- (if (is.matrix(tmp)) colSums(tmp, na.rm = TRUE) else tmp) >= 1
    Features <- features[featIds]
    models <- lapply(Features, function(gene) {
        df <- buildDataFrame(obj, gene = gene, pi = pi)
        if (is.null(df) || sum(!is.na(df$pi)) < 3) {
            return(NULL)
        }
        if (randomNested) {
            df <- nestRandom(df, randomVarsSplit)
        }
        if (MM) {
            mod <- try(lmerTest::lmer(Formula, data = df, na.action = na.omit, weights = if (noWeight) {
                NULL
            } else {
                df$weight
            }, contrasts = contrasts, control = Control), silent = TRUE)
        }
        # If no random variables or fit failed, go for linear model
        if (!MM || is(mod, "try-error")) {
            mod <- try(lm(formula(fixedPart), data = df, na.action = na.omit, weights = if (noWeight) {
                NULL
            } else {
                df$weight
            }, contrasts = contrasts), silent = TRUE)
        }
        # If still fails, drop fixed effects
        if (is(mod, "try-error")) {
            mod <- lm(formula("pi - 0.5 ~ 1"), data = df, na.action = na.omit, weights = if (noWeight) {
                NULL
            } else {
                df$weight
            })
        }
        return(mod)
    })
    names(models) <- Features
    results <- extractResults(models, hypFrame = obj$hypFrame, fixedVars)
    # Effect size, standard error, p-value and adjusted p-value per mixed effect
    if (!returnModels) {
        return(list(results = results))
    } else {
        return(list(results = results, models = models))
    }
}
#' Fit linear (mixed) models fitting for all probabilistic indices (PIs)
#' @description The PI is used as outcome variable in a linear (mixed) model,
#' with design variables as regressors.
#'
#' @param obj The result object
#' @param pis Optional, the pis required. Defaults to all pis in the object
#' @param verbose a boolean, should output be printed?
#' @param ... Passed onto fitLMMsSingle
#'
#' @return A list of fitted objects
#' @export
#'
#' @examples
#' example(addWeightFunction, "spatrans")
#' lmmModels <- fitLMMs(yangObj, fixedVars = "day", randomVars = "root")
fitLMMs <- function(obj, pis = obj$pis, verbose = TRUE, ...) {
    out <- lapply(pis, function(pi) {
        fitLMMsSingle(obj, pi = pi, verbose = verbose, ...)
    })
    names(out) <- pis
    return(out)
}
