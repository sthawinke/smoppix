#' Extract a data frame for a certain gene and PI from a fitted object
#' @description
#' Based on a fitted object, a dataframe with results for a certain feature and
#'  PI is built, e.g. in preparation for linear modelling.
#'
#' @param obj A results object. For distances to fixed objects, the result of a call to \link{estPis};
#' for nearest neighbour distances, the result of a call to \link{addWeightFunction}
#' @param gene A character string indicating the desired gene or gene pair (genes separated by double hyphens)
#' @param pi character string indicating the desired PI
#' @param piMat A data frame. Will be constructed if not provided, for internal use.
#' @param moransI A boolean, should Moran's I be calculated?
#' in the linear mixed model
#' @param weightMats List of weight matrices for Moran's I calculation.
#' @param pppDf Dataframe of point pattern-wise variables. It is precalculated
#' in fitLMMsSingle for speed, but will be newly constructed when not provided.
#' @param prepMat A preconstructed matrix with PI estimates to avoid looping
#' @inheritParams fitLMMs
#' @return A dataframe with estimated PIs and covariates
#' @export
#' @importFrom scam predict.scam
#' @seealso \link{addWeightFunction}, \link{buildMoransIDataFrame}
#' @examples
#' example(addWeightFunction, "smoppix")
#' dfUniNN <- buildDataFrame(yangObj, gene = "SmVND2", pi = "nn")
#' # Example analysis with linear mixed model
#' library(lmerTest)
#' mixedMod <- lmer(pi - 0.5 ~ day + (1 | root),
#'     weight = weight, data = dfUniNN,
#'     contrasts = list("day" = "contr.sum")
#' )
#' summary(mixedMod)
#' # Evidence for aggregation
buildDataFrame <- function(obj, gene, pi = c("nn", "nnPair", "edge", "centroid",
                               "nnCell", "nnPairCell"), piMat, moransI = FALSE,
                           numNNs = 8, weightMats, pppDf, prepMat) {
    pi <- match.arg(pi)
    if (missing(pppDf)) {
        pppDf <- as.data.frame(obj$hypFrame[, c("image", getPPPvars(obj)), drop = FALSE])
    }
    if (missing(piMat)) {
        stopifnot((lg <- length(gene)) %in% c(1, 2))
        foo <- checkPi(obj, pi)
        if (!(windowId <- (pi %in% c("edge", "centroid"))) && is.null(obj$Wfs[[pi]])) {
            stop("No weight function added yet, run addWeightFunction first!")
        }
        # Establish whether pi and gene match, and call separate functions
        cellId <- grepl("Cell", pi)
        if (cellId) {
            piSub <- sub("Cell", "", pi)
        }
        if (pairId <- grepl("Pair", pi)) {
            if (lg == 2) {
                gene <- paste(gene, collapse = "--")
            } else if (!grepl("--", gene)) {
                stop("Provide gene pair as character vector of length 2", " or separated by '--'.")
            }
        } else {
            if (lg != 1 || grepl("--", gene)) {
                stop("Provide a single gene for univariate PIs!")
            }
        }
        geneSplit <- sund(gene)
        misPrep <- missing(prepMat)
        # Check whether genes are present
        if (misPrep && any(id <- !(geneSplit %in% getFeatures(obj)))) {
            stop("PIs for features\n", geneSplit[id], "\nnot found in object")
        }
        piListNameInner <- if (windowId) {
            "windowDists"
        } else if (cellId) {
            "withinCellDists"
        } else {
            "pointDists"
        }
        if (cellId || windowId) {
            eventVars <- getEventVars(obj)
        }
        piDfs <- lapply(seq_len(nrow(obj$hypFrame)), function(n) {
            nList <- obj$hypFrame[n, , drop = TRUE]
            class(nList) <- "list" # Simplify things
            # Extract pi, and add counts and covariates
            df <- if (windowId) {
                piEst <- getGp(nList$pimRes[[piListNameInner]], gene)[[pi]]
                dfWin <- data.frame(pi = unlist(piEst), cell = rep(names(piEst),
                    times = vapply(piEst, FUN.VALUE = double(1), length)
                ))
                if (length(cellVars <- setdiff(eventVars, c("gene", "cell"))) &&
                    NROW(dfWin)) {
                    mat <- marks(nList$ppp, drop = FALSE)[match(dfWin$cell, Marks$cell), cellVars, drop = FALSE]
                    colnames(mat) <- cellVars
                    dfWin <- cbind(dfWin, mat)
                }
                return(dfWin)
            } else if (cellId) {
              idG <- if(pairId){
                marks(nList$ppp, drop = FALSE)$gene == geneSplit[1] | marks(nList$ppp, drop = FALSE)$gene == geneSplit[2]
              } else{
                marks(nList$ppp, drop = FALSE)$gene == geneSplit
              }
              Marks <- marks(nList$ppp[idG, ], drop = FALSE)
              Marks <- Marks[Marks$cell != "NA",]
              piEst <- if(misPrep){
                vapply(nList$pimRes[[piListNameInner]], FUN.VALUE = double(1), function(y) {
                  getGp(y[[piSub]], gene, notFoundReturn = NA)
                })
              } else {
                getGp(t(prepMat[names(nList$pimRes[[piListNameInner]]),, drop = FALSE]), gene, notFoundReturn = NA)
              }
                cellCovars <- Marks[, eventVars, drop = FALSE]
                cellCovars <- cellCovars[match(names(piEst), cellCovars$cell), setdiff(
                    colnames(cellCovars), "gene"), drop = FALSE]
                tabCell <- table(Marks$gene, Marks$cell)
                npVec <- vapply(geneSplit, FUN.VALUE = integer(ncol(tabCell)), function(x) {
                  getGp(tabCell, x, notFoundReturn = NA)
                })
                colnames(npVec) <- if (pairId) {
                    c("minP", "maxP")
                } else {
                    "NP"
                }
                npVecOut <- matrix(0, ncol = length(geneSplit), nrow = length(piEst), dimnames = list(names(piEst), colnames(npVec)))
                npVecOut[rownames(npVec),] <- npVec
                data.frame(pi = piEst, cellCovars, npVecOut)
            } else {
              piEst <- getGp(if(misPrep){nList$pimRes[[piListNameInner]][[pi]]} else {prepMat[n, ]},
                             gene, notFoundReturn = NA)
              npVec <- vapply(geneSplit, FUN.VALUE = integer(1), function(x) {
                getGp(nList$tabObs, x, notFoundReturn = NA)
              })
              if (pairId) {
                  cbind(pi = piEst, minP = min(npVec), maxP = max(npVec))
              } else {
                  cbind(pi = piEst, NP = npVec)
              }
            }
        })
        piMat <- if (!all((Times <- vapply(piDfs, FUN.VALUE = integer(1), NROW)) == 0)) {
          image <- rep(rownames(obj$hypFrame), times = Times)
          piDfsMat <- Reduce(piDfs, f = rbind)
          piMat <- data.frame(piDfsMat, pppDf[image, ,drop = FALSE])
          if (!windowId && !moransI) {
              # Add weights
              weight <- evalWeightFunction(obj$Wfs[[pi]], newdata = piMat[, if (grepl(
                  "Pair",
                  pi
              )) {
                  c("minP", "maxP")
              } else {
                  "NP"
              }, drop = FALSE])
              weight <- weight / sum(weight, na.rm = TRUE)
              piMat <- cbind(piMat, weight = weight)
          }
          rownames(piMat) <- NULL
          piMat
        }
      } else if (moransI) {
          if (missing(weightMats)) {
              weightMats <- lapply(getHypFrame(obj)$centroids, buildMoransIWeightMat,
                  numNNs = numNNs)
              names(weightMats) <- getHypFrame(obj)$image
          }
        piMat <- buildMoransIDataFrame(piMat = piMat, pi = pi, weightMats = weightMats)
      }
    return(piMat)
}
prepareMatrix <- function(obj, pi, features){
  if(windowId <- (pi %in% c("edge", "centroid"))){
    stop("prepareMatrix() not implemented yet for pi ", pi, " !")
  } else if (is.null(obj$Wfs[[pi]])) {
    stop("No weight function added yet, run addWeightFunction first!")
  }
  # Establish whether pi and gene match, and call separate functions
  if (cellId <- grepl("Cell", pi)) {
    piSub <- sub("Cell", "", pi)
  }
  piListNameInner <- if (windowId) {
    "windowDists"
  } else if (cellId) {
    "withinCellDists"
  } else {
    "pointDists"
  }
  if(missing(features)){
    features <- if(windowId){
      unique(unlist(lapply(obj$hypFrame$pimRes, function(x) names(x[[piListNameInner]]))))
    } else if(cellId){
      unique(unlist(lapply(obj$hypFrame$pimRes, function(x) lapply(x[[piListNameInner]], function(y) names(y[[piSub]])))))
    } else {
      unique(unlist(lapply(obj$hypFrame$pimRes, function(x){
        names(x[[piListNameInner]][[pi]])
      })))
    }
  }
  samples <- if(windowId){
    unique(unlist(lapply(obj$hypFrame$pimRes, function(x) lapply(x[[piListNameInner]], function(y) names(y[[pi]])))))
  } else if(cellId){
    unique(unlist(lapply(obj$hypFrame$pimRes, function(x) names(x[[piListNameInner]]))))
  } else {rownames(obj$hypFrame)}
  newMat <- matrix(NA, nrow = length(samples), ncol = length(features),
                  dimnames = list(samples, features))
  if(windowId){
    for(feat in features){
      featVec <- obj$hypFrame[[sam, "pimRes"]][[piListNameInner]][[feat]][[pi]]
      newMat[sam, features] <- featVec[features]
    }
  } else if(cellId){
    for(sam in rownames(obj$hypFrame)){
      samVec <- names(obj$hypFrame$pimRes[[sam]][[piListNameInner]])
      for(samIn in samVec){
        featVec <- obj$hypFrame$pimRes[[sam]][[piListNameInner]][[samIn]][[piSub]]
        newMat[samIn, features] <- featVec[features]
      }

    }
  } else {
    for(sam in samples){
      featVec <- obj$hypFrame[[sam, "pimRes"]][[piListNameInner]][[pi]]
      newMat[sam, features] <-featVec[features]
    }
  }
  return(newMat)
}
