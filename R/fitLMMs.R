#' Fit linear mixed nodels for all features of a pimRes object
#'
#' @inheritParams buildDfMM
#' @param fixedVars,randomVars Names of fixed and random variables
#' @param verbose A boolean, should the formula be printed?
#' @param returnModels a boolean: should the full models be returned?
#' Otherwise only summary statistics are returned
#' @param Formula A formula; if not supplied it will be constructed
#' from the fixed and random variables
#' @details Genes or gene pairs with insufficient observations will be silently
#' omitted
#'
#' @return A fitted linear mixed model of class 'lmerTest'
#' @export
#' @importFrom lmerTest lmer
#' @importFrom stats formula na.omit anova lm
#' @importFrom lme4 lmerControl .makeCC isSingular
#' @importFrom BiocParallel bplapply
#'
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPims(hypYang, pis = "nn",
#' features = attr(hypYang, "features")[seq_len(20)])
#' # First build the weighting function
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' lmmModels <- fitLMMs(yangObj, fixedVars = "day", randomVars = "root",
#' pi = "nn")
#' @seealso \link{buildDfMM},\link{getResults}
#' #FIX ME: write wrapper over all PIs!
fitLMMs <- function(resList, pi, fixedVars = NULL, randomVars = NULL,
                    verbose = TRUE, returnModels = FALSE, Formula = NULL) {
    pi <- match.arg(pi, choices = c(
        "nn", "allDist", "nnPair", "allDistPair",
        "edge", "midpoint", "nnCell", "allDistCell",
        "nnPairCell", "allDistPairCell"
    ))
    noWeight <- pi %in% c("edge", "midpoint")
    # For independent distances, no weights are needed
    designVars <- c(fixedVars, randomVars)
    stopifnot(all(designVars %in% c(resList$designVars,
                                    attr(resList$hypFrame, "cellVars"))))
    #Allow cell as design variable, both fixed and random
    if (is.null(Formula)) {
        fixedPart = paste(
            "pi - 0.5 ~ 1",
            paste(if (!is.null(fixedVars)) {
                    paste("+", paste(fixedVars, collapse = "+"))
                }))
        Formula <- formula(formChar <- paste(fixedPart, if(!is.null(randomVars)){
                    paste("+", paste0("(1|", randomVars, ")", collapse = "+"))
                }
            ))
    } else {
        formChar = paste(as.character(Formula), collapse = "")
    }
    MM <- any(grepl("\\|", formChar))
    if (verbose) {
        message("Fitted formula:\n", formChar)
    }
    Control <- lmerControl( ## convergence checking options
        check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL),
        check.conv.singular = .makeCC(action = "ignore",
                                      tol = formals(isSingular)$tol),
        check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6)
    )
    if (is.null(fixedVars)) {
        contrasts <- NULL
    } else {
        names(fixedVars) = fixedVars
        contrasts <- lapply(fixedVars, function(x) "named.contr.sum")
    }
    Features <- if (grepl("Pair", pi)) {
        feats = makePairs(attr(resList$hypFrame, "featuresEst"))
        featIds = colSums(vapply(feats, FUN.VALUE = double(nrow(resList$hypFrame)),
             function(gene) {
                 vapply(resList$hypFrame$tabObs, FUN.VALUE = double(1),
                        function(x) all(x[sund(gene)] >= 1))
             }), na.rm = TRUE) >= 1
        feats[featIds]
    } else {
        feats = attr(resList$hypFrame, "featuresEst")
        #First check if gene is present at all
        featIds = colSums(vapply(feats,
            FUN.VALUE = double(nrow(resList$hypFrame)), function(gene) {
                vapply(resList$hypFrame$tabObs, FUN.VALUE = double(1),
                       function(x) x[gene])
            }) >1, na.rm = TRUE) >= 1
        feats[featIds]
        #Leave out barely expressed genes
    }
    models <- bplapply(Features, function(gene) {
        df <- buildDfMM(resList, gene = gene, pi = pi)
        if (is.null(df) && sum(!is.na(df$pi)) < 3) {
            return(NULL)
        }
        if (MM) {
            mod <- try(lmerTest::lmer(Formula,
                data = df, na.action = na.omit,
                weights = if (noWeight) NULL else weight,
                contrasts = contrasts, control = Control
            ), silent = TRUE)
        }
        # If no random variables or fit failed, go for linear model
        if (!MM || is(mod, "try-error")) {
            mod <- try(lm(formula(fixedPart),
                data = df, na.action = na.omit,
                weights = if (noWeight) NULL else weight,
                contrasts = contrasts
            ), silent = TRUE)
        }
        #If still fails, drop fixed effects
        if (is(mod, "try-error")) {
            mod <- lm(formula("pi - 0.5 ~ 1"), data = df, na.action = na.omit,
                      weights = if (noWeight) NULL else weight
            )
        }
        return(mod)
    })
    names(models) = Features
    results <- extractResults(models, hypFrame = resList$hypFrame, fixedVars)
    # Effect size, standard error, p-value and adjusted p-value per mixed effect
    if (!returnModels) {
        return(list("results" = results, "resList" = resList))
    } else {
        return(list("results" = results, "resList" = resList,
                    "models" = models))
    }
}

