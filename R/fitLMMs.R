#' Fit linear (mixed) models for all probabilistic indices (PIs) and all genes
#' @description The PI is used as outcome variable in a linear (mixed) model,
#' with design variables as regressors. Separate models are fitted for every combination of gene
#' and PI. fitLMMsSingle() is the workhorse function for a single point pattern,
#'
#' @param obj The result object, from a call to \link{estPis} of \link{addWeightFunction}
#' @param pis Optional, the pis required. Defaults to all pis in the object
#' @param verbose a boolean, should output be printed?
#' @param fixedVars Names of fixed effects
#' @param randomVars Names of random variables
#' @param verbose A boolean, should the formula be printed?
#' @param returnModels a boolean: should the full models be returned?
#' Otherwise only summary statistics are returned
#' @param Formula A formula; if not supplied it will be constructed
#' from the fixed and random variables
#' @param randomNested A boolean, indicating if random effects are nested within
#' point patterns. See details.
#' @param features The features for which to fit linear mixed models.
#' Defaults to all features in the object
#' @param moranFormula Formula for Moran's I model fitting
#' @param addMoransI A boolean, include Moran's I of the cell-wise PIs in the calculation
#' @param numNNs An integer, the number of nearest neighbours in the
#' weight matrix for the calculation of the Moran's I statistic
#' @param ... Passed onto fitLMMsSingle
#' @details Genes or gene pairs with insufficient observations will be silently
#' omitted. When randomVars is provided as a vector, independent random
#' intercepts are fitted for them by default. Providing them separated by '\' or
#' ':' as in the lmer formulas is also allowed to reflect nesting structure, 
#' but the safest is to construct the formula yourself and pass it onto fitLMMs.
#'
#' It is by default assumed that random effects are nested within the point
#'  patterns. This means for instance that cells with the same name but from
#'  different point patterns are assigned to different random effects. Set
#'  'randomNested' to FALSE to override this behaviour.
#'
#'  The Moran's I statistic is used to test whether cell-wise PIs ("nnCell", "nnCellPair", "edge" and "centroid") 
#'  are spatially autocorrelated across the images. The numeric value of the PI is assigned to the 
#'  centroid location, and then Moran's I is calculated with a fixed number of numNNs nearest neighbours with equal weights.
#' @return For fitLMMs(), a list of fitted objects
#' @export
#' @seealso \link{buildMoransIDataFrame}, \link{buildDataFrame}
#' @order 1
#' @examples
#' example(addWeightFunction, "smoppix")
#' lmmModels <- fitLMMs(yangObj, fixedVars = "day", randomVars = "root")
#' res <- getResults(lmmModels, "nn", "Intercept") #Extract the results
#' head(res)
fitLMMs <- function(
    obj, pis = obj$pis, fixedVars = NULL, randomVars = NULL, verbose = TRUE,
    returnModels = FALSE, Formula = NULL, randomNested = TRUE, features = getEstFeatures(obj),
    moranFormula = NULL, addMoransI = FALSE, numNNs = 10, ...) {
  stopifnot(
    is.logical(returnModels), is.logical(randomNested), is.logical(addMoransI),
    is.numeric(numNNs)
  )
  if (addMoransI) {
    weightMats <- lapply(getHypFrame(obj)$centroids, buildMoransIWeightMat, numNNs = numNNs)
    names(weightMats) <- getHypFrame(obj)$image
  }
  out <- lapply(pis, function(pi) {
    fitLMMsSingle(obj, pi = pi, verbose = verbose, fixedVars = fixedVars, randomVars = randomVars,
        returnModels = returnModels, Formula = Formula, randomNested = randomNested,
        features = features, weightMats = weightMats, moranFormula = moranFormula,
        addMoransI = addMoransI, ...
    )
  })
  names(out) <- pis
  return(out)
}
#' @param weightMats A list of weight matrices for Moran's I
#' @return For fitLMMsSingle(), a list of test results, if requested also the linear models are returned
#' @importFrom lmerTest lmer
#' @importFrom stats formula anova
#' @importFrom lme4 lmerControl .makeCC isSingular
#' @importFrom methods is
#' @rdname fitLMMs
#' @order 2
fitLMMsSingle <- function(
    obj, pi, fixedVars, randomVars, verbose, returnModels,
    Formula, randomNested, features, addMoransI, weightMats, moranFormula) {
  pi <- match.arg(pi, choices = c(
    "nn", "nnPair", "edge", "centroid", "nnCell", "nnPairCell"))
  noWeight <- pi %in% c("edge", "centroid")
  if (noWeight) {
    randomVars <- setdiff(union(randomVars, "image/cell"), c("image", "cell"))
    # For edge and centroid, cell is nested within image
  }
  if ((cellId <- grepl("Cell", pi))) {
    randomVars <- union(randomVars, "image")
    # For cell, the image is always a design factor
  }
  if (pi %in% c("nn", "nnPair") && any(fixedVars %in% (evVars <- getEventVars(obj)))) {
    fixedVars <- setdiff(fixedVars, evVars)
    warning("Cell-wise variables cannot be incorporated into analysis with pi ",
            pi, ", so variables\n", paste(evVars, collapse = ", "), "\nwill be dropped.",
            immediate. = TRUE
    )
  }
  # For independent distances, no weights are needed
  randomVarsSplit <- if (!is.null(randomVars)) {
    grep("[[:punct:]]", value = TRUE, invert = TRUE, unique(unlist(lapply(c(
      "/",
      ":"
    ), function(Split) {
      lapply(randomVars, function(x) {
        strsplit(x, Split)[[1]]
      })
    }))))
  }
  designVars <- c(fixedVars, randomVarsSplit)
  if (any(id <- !(designVars %in% c(getDesignVars(obj), "image")))) {
    stop("Design variables ", paste(designVars[id], collapse = " "), " not found in object.")
  }
  Formula <- buildFormula(Formula, fixedVars, randomVars)
  MM <- any(grepl("\\|", formChar <- characterFormula(Formula)))
  if (verbose) {
    message("Fitted formula for pi ", pi, ":\n", formChar)
  }
  if (addMoransI) {
    moranFormula <- buildFormula(moranFormula,
                                 fixedVars = morFix <- intersect(
                                   fixedVars,
                                   getPPPvars(obj)
                                 ), randomVars = rvMoran <- intersect(randomVars, getPPPvars(obj)),
                                 outcome = "MoransI"
    )
    MMmoran <- as.logical(length(rvMoran))
    names(morFix) <- morFix
    contrastsMoran <- lapply(morFix, function(x) named.contr.sum)
    if (verbose) {
      message("Fitted formula for Moran's I for pi ", pi, ":\n", formCharMoran <- characterFormula(moranFormula))
    }
  }
  Control <- lmerControl(
    check.conv.grad = .makeCC("ignore", tol = 0.002, relTol = NULL),
    check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),
    check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)
  )
  # convergence checking options
  tmp <- if (grepl("Pair", pi)) {
    features <- makePairs(features)
    vapply(features, FUN.VALUE = double(nrow(obj$hypFrame)), function(gene) {
      vapply(obj$hypFrame$tabObs, FUN.VALUE = double(1), function(x) {
        all(x[sund(gene)] >= 1)
      })
    })
  } else {
    # First check if gene is present at all
    vapply(features, FUN.VALUE = double(nrow(obj$hypFrame)), function(gene) {
      vapply(obj$hypFrame$tabObs, FUN.VALUE = double(1), function(x) x[gene])
    })
  }
  featIds <- (if (is.matrix(tmp)) {
    colSums(tmp, na.rm = TRUE)
  } else {
    tmp
  }) >= 1
  Features <- features[featIds]
  pppDf <- centerNumeric(as.data.frame(obj$hypFrame[, c("image", getPPPvars(obj))]))
  if (is.null(fixedVars)) {
    contrasts <- NULL
  } else {
    discreteVars <- intersect(getDiscreteVars(obj), fixedVars)
    names(discreteVars) <- discreteVars
    contrasts <- lapply(discreteVars, function(x) named.contr.sum)
  }
  prepMat = prepareMatrix(obj, pi = pi)
  models <- loadBalanceBplapply(Features, function(gene) {
    df <- buildDataFrame(obj, gene = gene, pi = pi, pppDf = pppDf, prepMat = prepMat)
    out <- if (is.null(df) || sum(!is.na(df$pi)) < 3) {
      NULL
    } else {
      moranMod <- if (addMoransI) {
        finalDf <- buildDataFrame(piMat = df, gene = gene, pi = pi,
          moransI = TRUE, obj = obj, weightMats = weightMats)
        W <- 1 / finalDf$Variance
        fitPiModel(moranFormula, finalDf, contrastsMoran, Control,
                   MM = MMmoran, Weight = W / sum(W, na.rm = TRUE))
      } # Run this before nesting
      if (randomNested) {
        df <- nestRandom(df, randomVarsSplit, intersect(fixedVars, getPPPvars(obj)))
      }
      contrasts <- contrasts[!names(contrasts) %in% vapply(df, FUN.VALUE = TRUE, is.numeric)]
      piMod <- fitPiModel(Formula, df, contrasts,
                          Control, MM = MM, Weight = df$weight)
      list(piMod = piMod, moranMod = moranMod)}
    return(out)
  })
  names(models) <- Features
  mods <- c("piMod", if (addMoransI) "moranMod")
  names(mods) <- mods
  results <- lapply(mods, function(mm) {
    extractResults(models, hypFrame = obj$hypFrame, fixedVars, subSet = mm)
  })
  # Effect size, standard error, p-value and adjusted p-value per mixed effect
  return(list(results = results$piMod, resultsMoran = results$moranMod, models = if(returnModels) models))
}
