#' Fit linear mixed nodels for all features of a pimRes object
#'
#' @inheritParams buildDataFrame
#' @inheritParams fitLMMs
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
#' @importFrom stats formula anova
#' @importFrom lme4 lmerControl .makeCC isSingular
#' @importFrom methods is
#' @importFrom ape Moran.I
#' @seealso \link{buildDataFrame},\link{getResults}
fitLMMsSingle <- function(obj, pi, fixedVars, randomVars, verbose, returnModels,
                          Formula, randomNested, features, addMoransI, weightMats, moranFormula) {
    pi <- match.arg(pi, choices = c(
        "nn", "nnPair", "edge", "centroid",
        "nnCell", "nnPairCell"
    ))
    noWeight <- pi %in% c("edge", "centroid")
    if (noWeight) {
        randomVars <- setdiff(union(randomVars, "image/cell"), c("image", "cell"))
        #For edge and centroid, cell is nested within image
    }
    if ((cellId <- grepl("Cell", pi))) {
        randomVars <- union(randomVars, "image")
        #For cell, the image is always a design factor
    }
    if(pi %in% c("nn", "nnPair") && any(fixedVars %in% (evVars <- getEventVars(obj)))){
        fixedVars = setdiff(fixedVars, evVars)
        warning("Cell-wise variables cannot be incorporated into analysis with pi ",
                pi, ", so variables\n", paste(evVars, collapse = ", "),
                "\nwill be dropped.", immediate. = TRUE)
    }
    # For independent distances, no weights are needed
    randomVarsSplit <- if (!is.null(randomVars)) {
        grep("[[:punct:]]", value = TRUE, invert = TRUE,
        unique(unlist(lapply(c("/", ":"), function(Split) {
            lapply(randomVars, function(x){
                strsplit(x, Split)[[1]]
            })
        }))))
    }
    designVars <- c(fixedVars, randomVarsSplit)
    if (any(id <- !(designVars %in% c(getDesignVars(obj), "image")))) {
        stop("Design variables ", paste(designVars[id], collapse = " "),
             " not found in object.")
    }
    Formula = buildFormula(Formula, fixedVars, randomVars)
    MM <- any(grepl("\\|", formChar <- characterFormula(Formula)))
    if (verbose) {
        message("Fitted formula for pi ", pi, ":\n", formChar)
    }
    if(addMoransI){
        moranFormula = buildFormula(moranFormula, fixedVars = intersect(fixedVars, getPPPvars(obj)),
                                     randomVars = intersect(randomVars, getPPPvars(obj)), outcome = "MoransI")
        if (verbose) {
            message("Fitted formula for Moran's I for pi ", pi, ":\n", formCharMoran <- characterFormula(moranFormula))
        }
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
            vapply(obj$hypFrame[, "tabObs", drop = TRUE], FUN.VALUE = double(1), function(x) all(x[sund(gene)] >= 1))
        })
    } else {
        # First check if gene is present at all
        vapply(features, FUN.VALUE = double(nrow(obj$hypFrame)), function(gene) {
            vapply(obj$hypFrame[, "tabObs", drop = TRUE], FUN.VALUE = double(1), function(x) x[gene])
        })
    }
    featIds <- (if (is.matrix(tmp)) colSums(tmp, na.rm = TRUE) else tmp) >= 1
    Features <- features[featIds]
    models <- loadBalanceBplapply(Features, function(gene) {
        df <- buildDataFrame(obj, gene = gene, pi = pi)
        if (is.null(df) || sum(!is.na(df$pi)) < 3) {
            return(NULL)
        }
        moranMod = if(addMoransI){
            dfMoran <- if(pi %in% c("edge", "centroid")){
                aggregate(pi ~ cell + image, df, FUN = mean, na.rm = TRUE)
            } else df
            moransIs = vapply(unIm <- unique(dfMoran$image), FUN.VALUE = double(1), function(im){
                subDf = dfMoran[dfMoran$image == im, ]
                moranObj = Moran.I(subDf$pi, weight = weightMats[[im]][subDf$cell, subDf$cell])
                (moranObj$observed-moranObj$expected)/moranObj$sd
            })
            finalDf = cbind("MoransI" = moransIs, df[match(unIm, df$image),
                                                          setdiff(colnames(df), "weight")])
            fitPiModel(moranFormula, finalDf, contrasts, Control, MM = MM)
        } #Run this before nesting
        if (randomNested) {
            df <- nestRandom(df, randomVarsSplit, intersect(fixedVars, getPPPvars(obj)))
        }
        piMod = fitPiModel(Formula, df, contrasts, Control, MM = MM)
        return(list("piMod" = piMod, "moranMod" = moranMod))
    })
    names(models) <- Features
    mods = c("piMod", if(addMoransI) "moranMod");names(mods) = mods
    results<- lapply(mods, function(mm){
        extractResults(models, hypFrame = obj$hypFrame, fixedVars, subSet = mm)
        })
    # Effect size, standard error, p-value and adjusted p-value per mixed effect
    if (!returnModels) {
        return(list(results = results$piMod, resultsMoran = results$moranMod))
    } else {
        return(list(results = results$piMod, resultsMoran = results$moranMod, models = models))
    }
}
#' Fit linear (mixed) models fitting for all probabilistic indices (PIs)
#' @description The PI is used as outcome variable in a linear (mixed) model,
#' with design variables as regressors.
#'
#' @param obj The result object
#' @param pis Optional, the pis required. Defaults to all pis in the object
#' @param verbose a boolean, should output be printed?
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
#' @param ... Passed onto fitLMMsSingle
#'
#' @return A list of fitted objects
#' @export
#'
#' @examples
#' example(addWeightFunction, "spatrans")
#' lmmModels <- fitLMMs(yangObj, fixedVars = "day", randomVars = "root")
fitLMMs <- function(obj, pis = obj$pis, fixedVars = NULL, randomVars = NULL,
                    verbose = TRUE, returnModels = FALSE, Formula = NULL,
                    randomNested = TRUE, features = getFeatures(obj), moranFormula = NULL,
                    addMoransI = any(pis %in% c("edge", "centroid", "nnCell", "nnPairCell")), ...){
    if(addMoransI){
        weightMats = lapply(getHypFrame(obj)$centroids, buildWeightMat)
    }
    out <- lapply(pis, function(pi) {
        fitLMMsSingle(obj, pi = pi, verbose = verbose, fixedVars = fixedVars,
                      randomVars = randomVars, returnModels = returnModels,
                      Formula = Formula, randomNested = randomNested,
                      features = features, weightMats = weightMats, moranFormula = moranFormula,
                      addMoransI = pi %in% c("edge", "centroid", "nnCell", "nnPairCell"), ...)
    })
    names(out) <- pis
    return(out)
}
#' @importFrom stats dist
buildWeightMat = function(coordList){
    coordMat = t(vapply(coordList, FUN.VALUE = double(2), getCoordsMat))
    as.matrix(1/dist(coordMat))
}

