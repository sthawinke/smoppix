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
#' @return For fitLMMs(), a list of fitted objects
#' @export
#' @seealso \link{buildDataFrame}
#' @order 1
#' @examples
#' example(addWeightFunction, "smoppix")
#' lmmModels <- fitLMMs(yangObj, fixedVars = "day", randomVars = "root")
#' res <- getResults(lmmModels, "nn", "Intercept") #Extract the results
#' head(res)
fitLMMs <- function(
    obj, pis = obj$pis, fixedVars = NULL, randomVars = NULL, verbose = TRUE,
    returnModels = FALSE, Formula = NULL, randomNested = TRUE, features = getEstFeatures(obj), ...) {
  stopifnot(is.logical(returnModels), is.logical(randomNested))
  out <- lapply(pis, function(pi) {
    fitLMMsSingle(obj, pi = pi, verbose = verbose, fixedVars = fixedVars, randomVars = randomVars,
        returnModels = returnModels, Formula = Formula, randomNested = randomNested,
        features = features, ...
    )
  })
  names(out) <- pis
  return(out)
}
#' @return For fitLMMsSingle(), a list of test results, if requested also the linear models are returned
#' @importFrom lmerTest lmer
#' @importFrom stats formula terms model.matrix
#' @importFrom lme4 lmerControl .makeCC isSingular lFormula mkReTrms findbars
#' @importFrom methods is
#' @rdname fitLMMs
#' @order 2
fitLMMsSingle <- function(
    obj, pi, fixedVars, randomVars, verbose, returnModels,
    Formula, randomNested, features) {
  pi <- match.arg(pi, choices = c(
    "nn", "nnPair", "edge", "centroid", "nnCell", "nnPairCell"))
  foo <- checkPi(obj, pi)
  if (windowId <- (pi %in% c("edge", "centroid"))) {
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
  if(!windowId){
    prepMatorList <- prepareMatrixOrList(obj, pi = pi, features = Features)
  }
  if(windowId){
    models <- loadBalanceBplapply(Features, function(gene) {
      df <- buildDataFrame(obj, gene = gene, pi = pi, pppDf = pppDf)
      out <- if (is.null(df) || sum(!is.na(df$pi)) < 3) {
        NULL
      } else {
        if (randomNested) {
          df <- nestRandom(df, randomVarsSplit, intersect(fixedVars, getPPPvars(obj)))
        }
        contrasts <- contrasts[!names(contrasts) %in% vapply(df, FUN.VALUE = TRUE, is.numeric)]
        fitPiModel(Formula, df, contrasts,
                            Control, MM = MM, Weight = df$weight)
      }
      return(out)
    })
  } else {
    if(cellId){
      prepTableOrList <- lapply(obj$hypFrame$ppp, function(x){
        tab <- table(marks(x)[, c("gene", "cell")])
        class(tab) = "matrix"
        t(tab[, colnames(tab)!= "NA"])
      })
      prepCells <- lapply(seq_along(prepTableOrList), function(n){
        eventMarks <- marks(obj$hypFrame$ppp[[n]], drop = FALSE)[, getEventVars(obj), drop = FALSE]
        eventMarks[match(colnames(prepTableOrList[[n]]), eventMarks$cell),]
      })
    } else {
      prepTableOrList = {
        singleFeats = getFeatures(obj)
        emptyTab = matrix(NA, nrow = nrow(pppDf), ncol = length(singleFeats), 
                          dimnames = list(rownames(pppDf), singleFeats))
        for(i in rownames(pppDf)){
          emptyTab[i,names(getHypFrame(obj)[[i,"tabObs"]])] = getHypFrame(obj)[[i,"tabObs"]]
        }
        emptyTab
      }
    }
    #Prepare generic dataframe with fixed and random effects, later 
    # swap in weights and outcome per feature
    baseDf = if(cellId){
      Reduce(f = rbind, lapply(seq_along(prepCells), function(n) {
        cbind(prepCells[[n]][, setdiff(colnames(prepCells[[n]]), "gene"), drop = FALSE],
              pppDf[n,setdiff(colnames(pppDf), "cell"), drop = FALSE])
      }))
    } else {
      pppDf
    }
    contrasts <- contrasts[!names(contrasts) %in% vapply(baseDf, FUN.VALUE = TRUE, is.numeric)]
    ff <- lFormula(Formula, data = data.frame("pi" = 0.5, baseDf), 
                   contrasts = contrasts, na.action = na.omit)
    Terms = terms(Formula)
    Attr = attr(ff$X, "assign")
    modMat = model.matrix(formula(paste("~", paste(collapse = "+", fixedVars))),
                          pppDf, contrasts.arg = contrasts) #Fixed effects model matrix
    if (randomNested) {
      ff$fr <- nestRandom(ff$fr, randomVarsSplit, intersect(fixedVars, getPPPvars(obj)))
      ff$reTrms <- mkReTrms(findbars(Formula[[length(Formula)]]), ff$fr)
    }
    #For cell wise properties, also include prepTabs and prepCells
    models <- loadBalanceBplapply(Features, function(gene) {
      mat = getPiAndWeights(obj, gene = gene, pi = pi, prepMat = prepMatorList,
                            prepTab = prepTableOrList)
      out <- if (is.null(mat) || sum(id <- !is.na(mat[, "pi"])) < 3) {
        NULL
      } else {
        ff$fr = ff$fr[id,, drop = FALSE];ff$X = ff$X[id,, drop = FALSE]
        attr(ff$X, "assign") = Attr
        ff$reTrms$Zt = ff$reTrms$Zt[, id, drop = FALSE]
        fitSingleLmmModel(ff = ff, y = mat[id, "pi"], Terms = Terms, modMat = modMat[id,],
                    weights = mat[id, "weights"], Control = Control)
        }
      return(out)
    })
  }
  names(models) <- Features
  results <- extractResults(models, hypFrame = obj$hypFrame, fixedVars)
  # Effect size, standard error, p-value and adjusted p-value per mixed effect
  return(list(results = results, models = if(returnModels) models))
}
#' Take an existing frame, add outcome and weight and fit lmer model
#'
#' @param ff The prepared frame
#' @param y outcome vector
#' @param weights weights vector
#'
#' @returns A fitted lmer model
#' @importFrom lme4 mkLmerDevfun optimizeLmer mkMerMod
#' @importFrom lmerTest as_lmerModLmerTest
fitSingleLmmModel <- function(ff, y, Control, Terms, modMat, weights = NULL) {
  fr <- ff$fr                    # this is a data.frame (model frame)
  ## Use model-frame column names used by stats::model.frame
  fr$`(weights)` = weights
  fr[["pi - 0.5"]] <- y - 0.5               # replace response
  id = !is.na(y)
  
  mod = try({
    devfun <- mkLmerDevfun(fr, ff$X, ff$reTrms, control = Control)
    opt <- optimizeLmer(devfun, control = Control)
    out <- mkMerMod(rho = environment(devfun), opt = opt, 
                    reTrms = ff$reTrms, fr = fr)
    out <- lmerTest:::as_lmerModLT(out, devfun = devfun)
  }, silent = TRUE)
  # Switch to fixed effects model when fit failed
  if(inherits(mod, "try-error")){
    lm_from_wfit(lm.wfit(y = y, x = modMat, w = weights), y = y, Terms = Terms)
  }
  return(mod)
}
lm_from_wfit <- function(obj, y, Terms, weights = NULL) {
  # Create a minimal lm object
  obj <- c(obj, list(y = y, terms = Terms))
  class(obj) <- "lm"
  return(obj)
}