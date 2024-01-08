#' Add a variance weighting function
#' @description Build a weighting function based on the data, and attach it to the object
#'
#' @param resList A results list, from a call to estPims()
#' @param designVars A character vector containing all design factors
#' (both fixed and random), that are also present as variables in hypFrame
#' @param ... Additional arguments passed on to the scam::scam function
#' @param maxObs,maxFeatures The maximum number of observations respectively
#' features for fitting the weighting function. See details
#' @param pis The PIs for which weighting functions are constructed
#' @param lowestLevelVar The design variable at the lowest level of nesting,
#' often separating technical replicates. The conditional variance is calculated
#' within the groups of PIs defined by this variable
#' @param minNumVar The minimum number of observations needed to calculate a
#' variance
#' @details For computational and memory reasons,
#'  for large datasets the trend fitting
#' is restricted to a subset of the data through the maxObs and maxFeatures
#' parameters.
#' Provide either 'designVars' or 'lowestLevelVar'. In the latter case, the
#' design variables are set to all imageVars in the hypFrame object except
#' lowestLevelVar. When the PI is calculated on the cell level,
#' this the lowest nesting level so inputs will be ignored for these PIs.
#' The registered parallel backend will be used for fitting the trends of the
#' different PIs.
#'
#' @return A fitted scam object, which can be used to
#' @details The scam functions fits a decreasing spline of the variance as a
#'  function of the number of observations
#' @importFrom scam scam
#' @importFrom stats formula
#' @importFrom Rfast rowAll
#' @export
#' @seealso \link{buildDataFrame}
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c('x', 'y'),
#'     imageVars = c('day', 'root', 'section')
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPims(hypYang[c(seq_len(5), seq(25, 29)),], nPointsAll = 1e3,
#'     pis = c('nn', 'nnPair'), features = getFeatures(hypYang)[1:10]
#' )
#' # Build the weighting function for all PIs present
#' yangObj <- addWeightFunction(yangPims, designVars = c('day', 'root'))
#' # Alternative formulation with 'lowestLevelVar'
#' yangObj2 <- addWeightFunction(yangPims, lowestLevelVar = 'section',
#'     pi = 'nn'
#' )
addWeightFunction <- function(resList, pis = resList$pis, designVars,
                              lowestLevelVar, maxObs = 1e+05, maxFeatures = 500,
                              minNumVar = 3, ...) {
    if (is.null(resList$pis)) {
        stop("No pims found in the hyperframe.",
             "First estimate them using the estPims() function.")
    }
    if (all(pis %in% c("edge", "midpoint"))) {
        stop("Calculating weight matrices for distances to fixed points is ",
        "unnecessary as they are independent.
             Simply proceed with fitting the model on the",
            " individual evaluations of the B-function.")
    }
    if(!is.null(resList$Wfs)){
        warning("Overwriting pre-existing weight function!")
    }
    pis <- match.arg(pis, choices = c("nn", "nnPair", "nnCell", "nnPairCell"), several.ok = TRUE)
    isNested <- length(getPPPvars(resList)) >= 1  #Is there a nested structure in the design
    allCell <- all(grepl("Cell", pis))
    if (isNested) {
        designVars <- constructDesignVars(designVars, lowestLevelVar, allCell, resList = resList)
    } else {
        # Only check design if there is nesting
        designVec <- integer(nrow(resList$hypFrame))
        designVars <- NULL
    }
    Wfs <- lapply(pis, function(pi) {
        windowId <- pi %in% c("edge", "midpoint")
        cellId <- any(grepl("Cell", pi))
        features <- getFeatures(resList)
        if (pairId <- grepl("Pair", pi)) {
            features <- makePairs(features)
        }
        if (cellId || !isNested) {
            ordDesign <- seq_len(nrow(resList$hypFrame))
        } else {
            designVec <- apply(as.data.frame(resList$hypFrame[, designVars, drop = FALSE]), 1, paste, collapse = "_")
            ordDesign <- order(designVec)  # Ensure correct ordering for tapply
        }
        varEls <- lapply(features, function(gene) {
            geneSplit <- if (pairId) {sund(gene)} else {gene}
            if (cellId) {
                piSub <- sub("Cell", "", pi)
                piList <- lapply(resList$hypFrame$pimRes, function(x) {
                  lapply(x[["withinCellDists"]], function(y) {
                    if (gene %in% names(y[[piSub]]))
                      y[[piSub]][[gene]]
                  })
                })
                # If cellId, there is no tapply, cells are the lowest level anyway
                tmp <- lapply(ordDesign, function(x) {
                  tab <- table(marks(resList$hypFrame$ppp[[x]])$cell,
                               marks(resList$hypFrame$ppp[[x]])$gene)
                  if (!all(geneSplit %in% colnames(tab))) {
                    return(NULL)
                  }
                  lenOut <- sum(id <- rowAll((CellGene <- tab[,
                                    geneSplit, drop = FALSE]) > (1 - pairId)))
                  if (!lenOut) {
                    return(NULL)
                  } else {
                    xx <- unlist(piList[[x]])
                    deps <- if (sum(!is.na(xx)) >= minNumVar) {
                      (xx - mean(xx, na.rm = TRUE))^2
                    } else {
                      rep_len(NA, lenOut)
                    }
                    cbind(quadDeps = matrix(deps, nrow = lenOut),
                          rowSort(CellGene[id, , drop = FALSE]))
                  }
                })
                if (all(vapply(tmp, FUN.VALUE = logical(1), is.null))) {
                  return(NULL)
                } else {
                  out <- t(Reduce(tmp, f = rbind))
                }
            } else {
                # Points with nn
                piList <- vapply(resList$hypFrame$pimRes[ordDesign],
                                 FUN.VALUE = double(1), function(x) {
                    if (is.null(baa <- getGp(x[["pointDists"]][[pi]], gene)))
                        NA else baa
                })
                quadDeps <- unlist(tapply(piList, designVec, function(x) {
                  if (sum(!is.na(x)) >= minNumVar) {
                    (x - mean(x, na.rm = TRUE))^2
                  } else {
                    rep_len(NA, length(x))
                  }
                }))  # The quadratic departures from the conditional mean
                tabEntries <- vapply(resList$hypFrame$tabObs[ordDesign],
                                     FUN.VALUE = double(if(pairId) 2 else 1),
                                     function(x) {
                    if (all(geneSplit %in% names(x))){
                      sort(x[geneSplit]) #Number of NN distances
                  } else {
                    rep(NA, if(pairId) 2 else 1)
                  }
                })
                out <- rbind("quadDeps" = quadDeps, tabEntries)
            }
            rownames(out)[-1] <- if(pairId) c("minP", "maxP") else "NP"
            out
        })
        varElMat <- matrix(unlist(varEls), ncol = if(pairId) 3 else 2, byrow = TRUE,
                           dimnames = list(NULL, c("quadDeps", if(pairId) c("minP", "maxP") else "NP")))
        # Faster than cbind
        varElMat <- varElMat[!is.na(varElMat[, "quadDeps"]) &
                                 varElMat[, "quadDeps"] != 0, ]
        if (nrow(varElMat) > maxObs) {
            varElMat <- varElMat[sample(nrow(varElMat), maxObs), ]
        }
        scamForm <- formula(paste("log(quadDeps) ~", if (pairId) {
            "s(log(minP), bs = 'mpd') + s(log(maxP), bs = 'mpd')"
        } else {"s(log(NP), bs = 'mpd')"}))
        scamMod <- scam(scamForm, data = data.frame(varElMat), ...)
        attr(scamMod, "pis") <- pi
        return(scamMod)
    })
    names(Wfs) <- pis
    return(c(resList[setdiff(names(resList), c("Wfs", "designVars"))],
             list(Wfs = Wfs, designVars = designVars)))
    #Overwrite preexisting Wfs
}
