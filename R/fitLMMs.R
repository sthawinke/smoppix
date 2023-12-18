#' Fit linear mixed nodels for all features of a pimRes object
#'
#' @inheritParams buildDfMM
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
#' intercepts are fitted for them by default. Providing them separated by "\" or
#' ":" as in the lmer formulas is also allowed to reflect nesting structure.
#'
#' It is by default assumed that random effects are nested within the point
#'  patterns. This means for instance that cells with the same name but from
#'  different point patterns are assigned to different random effects. Set
#'  "randomNested" to FALSE to override this behaviour.
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
#' yangPims <- estPims(hypYang[c(seq_len(5), seq(25, 29)),], pis = "nn")
#' # First build the weighting function
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' lmmModels <- fitLMMs(yangObj, fixedVars = "day", randomVars = "root",
#'     pi = "nn")
#' @seealso \link{buildDfMM},\link{getResults}
fitLMMs <- function(obj, pi, fixedVars = NULL, randomVars = NULL,
                    verbose = TRUE, returnModels = FALSE, Formula = NULL,
                    randomNested = TRUE, features = getFeatures(obj)) {
    pi <- match.arg(pi, choices = c(
        "nn", "allDist", "nnPair", "allDistPair",
        "edge", "midpoint", "nnCell", "allDistCell",
        "nnPairCell", "allDistPairCell"
    ))
    noWeight <- pi %in% c("edge", "midpoint")
    # For independent distances, no weights are needed
    randomVarsSplit = if(!is.null(randomVars))
        unlist(lapply(c("/", ":"), function(Split){
            tmp = strsplit(randomVars, Split)
            if(length(tmp)>2)
                tmp
        }))
    designVars <- c(fixedVars, randomVarsSplit)
    if(any(id <- !(designVars %in% getDesignVars(obj)))){
        stop("Design variables", designVars[id], "not found in object.")
    }
    # Allow cell as design variable, both fixed and random
    if (is.null(Formula)) {
        fixedPart <- paste(
            "pi - 0.5 ~ 1",
            paste(if (!is.null(fixedVars)) {
                paste("+", paste(fixedVars, collapse = " + "))
            })
        )
        Formula <- formula(formChar <- paste(fixedPart, if (!is.null(randomVars)) {
            paste("+", paste0("(1|", randomVars, ")", collapse = " + "))
        }))
    } else {
        formChar <- paste(as.character(Formula), collapse = "")
    }
    MM <- any(grepl("\\|", formChar))
    if (verbose) {
        message("Fitted formula:\n", formChar)
    }
    Control <- lmerControl( ## convergence checking options
        check.conv.grad = .makeCC("ignore", tol = 2e-3, relTol = NULL),
        check.conv.singular = .makeCC(action = "ignore",
            tol = formals(isSingular)$tol),
        check.conv.hess = .makeCC(action = "ignore", tol = 1e-6)
    )
    if (is.null(fixedVars)) {
        contrasts <- NULL
    } else {
        names(fixedVars) <- fixedVars
        contrasts <- lapply(fixedVars, function(x) named.contr.sum)
    }
    Features <- if (grepl("Pair", pi)) {
        feats <- makePairs(features)
        featIds <- colSums(vapply(feats,
            FUN.VALUE = double(nrow(obj$hypFrame)),
            function(gene) {
                vapply(obj$hypFrame$tabObs,
                    FUN.VALUE = double(1),
                    function(x) all(x[sund(gene)] >= 1)
                )
            }
        ), na.rm = TRUE) >= 1
        feats[featIds]
    } else {
        # First check if gene is present at all
        featIds <- colSums(vapply(features,
            FUN.VALUE = double(nrow(obj$hypFrame)), function(gene) {
                vapply(obj$hypFrame$tabObs,
                    FUN.VALUE = double(1),
                    function(x) x[gene]
                )
            }
        ) > 1, na.rm = TRUE) >= 1
        features[featIds]
        # Leave out barely expressed genes
    }
    models <- bplapply(Features, function(gene) {
        df <- buildDfMM(obj, gene = gene, pi = pi)
        if(randomNested){
            df = nestRandom(df, randomVarsSplit)
        }
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
        # If still fails, drop fixed effects
        if (is(mod, "try-error")) {
            mod <- lm(formula("pi - 0.5 ~ 1"),
                data = df, na.action = na.omit,
                weights = if (noWeight) NULL else weight
            )
        }
        return(mod)
    })
    names(models) <- Features
    results <- extractResults(models, hypFrame = obj$hypFrame, fixedVars)
    # Effect size, standard error, p-value and adjusted p-value per mixed effect
    if (!returnModels) {
        return(list("results" = results, "obj" = obj))
    } else {
        return(list(
            "results" = results, "obj" = obj,
            "models" = models
        ))
    }
}
#' A wrapper function for linear mixed model fitting for all PIs
#'
#' @param obj The result object
#' @param pis Optional, the pis required. Defaults to all pis in the object
#' @param ... Passed onto fitLMMs
#'
#' @return A list of fitted objects
#' @export
#'
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang, coordVars = c("x", "y"),
#' imageVars = c("day", "root", "section"))
#' yangPims <- estPims(hypYang, features = getFeatures(hypYang)[seq_len(10)],
#' pis = c("nn", "nnPair"))
#' allModels = fitLMMsWrapper(yangPims)
fitLMMsAll = function(obj, pis = obj$pis, ...){
    out = lapply(pis, function(pi){
        fitLMMs(obj, pi = pi, verbose = pi==pis[1],...)
    })
    names(out) = pis
    return(out)
}
