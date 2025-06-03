#' @description addWeightFunction() adds a weighting function based on the data to the object by modeling variance as a
#' non-increasing spline as a function of the number of events.
#'
#' @param resList A results list, from a call to estPis().
#' @param designVars A character vector containing all design factors
#' (both fixed and random), that are also present as variables in hypFrame.
#' @param ... Additional arguments passed on to the \link[scam]{scam} function, fitting the spline
#' @param maxObs,maxFeatures The maximum number of observations respectively
#' features for fitting the weighting function. See details.
#' @param pis The PIs to be estimated or for which weighting functions is to be added
#' @param lowestLevelVar The design variable at the lowest level of nesting,
#' often separating technical replicates. The conditional variance is calculated
#' within the groups of PIs defined by this variable.
#' @param minNumVar The minimum number of observations needed to calculate a
#' variance. Groups with fewer replicates are ignored.
#' @details Provide either 'designVars' or 'lowestLevelVar'. The 'designVars' are 
#' usually the same as the regressors in the linear model. In case 'lowestLevelVar' is provided, the
#' design variables are set to all imageVars in the hypFrame object except
#' lowestLevelVar. When the PI is calculated on the cell level ("nnCell" or "nnPairCell"),
#' the cell is always the lowest nesting level, and inputs to 'designVars' or 'lowestLevelVar' will be ignored for these PIs.
#' The registered parallel backend will be used for fitting the trends of the
#' different PIs.
#' For computational and memory reasons, for large datasets the trend fitting
#' is restricted to a random subset of the data through the maxObs and
#' maxFeatures parameters.
#'
#' @return For addWeightFunction(), the input object 'resList' with a slot 'Wfs' added containing the weighting
#' functions.
#' @importFrom scam scam
#' @importFrom stats formula
#' @export
#' @seealso \link{buildDataFrame}, \link{estPis}
#' @rdname estPis
#' @order 3
addWeightFunction <- function(resList, pis = resList$pis, designVars, lowestLevelVar,
    maxObs = 1e+05, maxFeatures = 1000, minNumVar = 3, ...) {
    if (is.null(resList$pis)) {
        stop("No pims found in the hyperframe.", " First estimate them using the estPis() function.")
    }
    if (all(pis %in% c("edge", "centroid"))) {
        stop(
            "Calculating weight matrices for distances to fixed points is ", "unnecessary as they are independent.
             Simply proceed with fitting the model on the",
            " individual evaluations of the B-function."
        )
    }
    if (!all(vapply(resList$Wfs[pis], FUN.VALUE = TRUE, is.null))) {
        warning("Overwriting pre-existing weight function!")
    }
    pis <- match.arg(pis, choices = c("nn", "nnPair", "nnCell", "nnPairCell"), several.ok = TRUE)
    isNested <- !(missing(designVars) && missing(lowestLevelVar)) # Is there a nested structure in the design
    allCell <- all(grepl("Cell", pis))
    if (isNested) {
        designVars <- constructDesignVars(designVars, lowestLevelVar, allCell, resList = resList)
    } else {
        # Only check design if there is nesting
        designVec <- integer(nrow(resList$hypFrame))
        designVars <- NULL
    }
    Wfs <- lapply(pis, function(pi) {
        cellId <- grepl("Cell", pi)
        features <- getEstFeatures(resList)
        if (pairId <- grepl("Pair", pi)) {
            features <- makePairs(features)
        }
        if (length(features) > maxFeatures) {
            features <- sample(features, maxFeatures)
        }
        if (cellId || !isNested) {
            ordDesign <- seq_len(nrow(resList$hypFrame))
        } else {
            designVec <- apply(as.data.frame(resList$hypFrame[, designVars, drop = FALSE]),
                1, paste, collapse = "_")
            ordDesign <- order(designVec) # Ensure correct ordering for tapply
        }
        varEls <- lapply(features, function(gene) {
            geneSplit <- if (pairId) {
                sund(gene)
            } else {
                gene
            }
            if (cellId) {
                piSub <- sub("Cell", "", pi)
                piList <- lapply(resList$hypFrame$pimRes, function(x) {
                    lapply(x[["withinCellDists"]], function(y) {
                        if (!is.null(vec <- getGp(y[[piSub]], gene))) {
                            vec
                        } else {
                            NULL
                        }
                    })
                })
             # If cellId, there is no tapply, cells are the lowest level anyway
                tmp <- lapply(ordDesign, function(x) {
                    tab <- table(marks(resList$hypFrame$ppp[[x]])$cell, marks(resList$hypFrame$ppp[[x]])$gene)
                    tab <- tab[setdiff(rownames(tab), "NA"), ]
                    out <- if (all(geneSplit %in% colnames(tab))) {
                      lenOut <- sum(id <- apply((CellGene <- tab[, geneSplit, drop = FALSE]) >
                          (1 - pairId), 1, all))
                      if (lenOut) {
                          xx <- unlist(piList[[x]])
                          deps <- if (sum(!is.na(xx)) >= minNumVar) {
                              (xx - mean(xx, na.rm = TRUE))^2
                          } else {
                              rep_len(NA, lenOut)
                          }
                          cbind(quadDeps = matrix(deps, nrow = lenOut), if (pairId) {
                              t(apply(CellGene[id, , drop = FALSE], 1, sort))
                          } else {
                              CellGene[id, , drop = FALSE]
                          })
                      }
                    }
                    return(out)
                })
                if (all(vapply(tmp, FUN.VALUE = logical(1), is.null))) {
                    return(NULL)
                } else {
                    out <- t(Reduce(tmp, f = rbind))
                }
            } else {
                # Points with nn
                piList <- vapply(resList$hypFrame$pimRes[ordDesign],
                    FUN.VALUE = double(1),
                    function(x) {
                        if (is.null(baa <- getGp(x[["pointDists"]][[pi]], gene))) {
                            NA
                        } else {
                            baa
                        }
                    }
                )
                quadDeps <- unlist(tapply(piList, designVec, function(x) {
                    if (sum(!is.na(x)) >= minNumVar) {
                        (x - mean(x, na.rm = TRUE))^2
                    } else {
                        rep_len(NA, length(x))
                    }
                })) # The quadratic departures from the conditional mean
                tabEntries <- vapply(resList$hypFrame$tabObs[ordDesign], FUN.VALUE = double(if (pairId) {
                    2
                } else {
                    1
                }), function(x) {
                    if (all(geneSplit %in% names(x))) {
                        sort(x[geneSplit]) # Number of NN distances
                    } else {
                        rep(NA, if (pairId) 2 else 1)
                    }
                })
                out <- rbind(quadDeps = quadDeps, tabEntries)
            }
            rownames(out)[-1] <- if (pairId) {
                c("minP", "maxP")
            } else {
                "NP"
            }
            out
        })
        varElMat <- matrix(unlist(varEls), ncol = if (pairId) {
            3
        } else {
            2
        }, byrow = TRUE, dimnames = list(NULL, c("quadDeps", if (pairId) {
            c("minP", "maxP")
        } else {
            "NP"
        })))
        # Build matrix with variance entries and number of events
        varElMat <- varElMat[!is.na(varElMat[, "quadDeps"]) &
                                 varElMat[, "quadDeps"] != 0, ]
        if (nrow(varElMat) > maxObs) {
            varElMat <- varElMat[sample(nrow(varElMat), maxObs), ]
        }
        scamForm <- formula(paste("log(quadDeps) ~", if (pairId) {
            "s(log(minP), bs = 'mpd') + s(log(maxP), bs = 'mpd')"
        } else {
            "s(log(NP), bs = 'mpd')"
        }))
        # Fit a spline
        scamMod <- scam(scamForm, data = data.frame(varElMat), ...)
        attr(scamMod, "pis") <- pi
        return(scamMod)
    })
    names(Wfs) <- pis
    return(c(resList[setdiff(names(resList), c("Wfs", "designVars"))], list(
        Wfs = Wfs,
        designVars = designVars
    )))
}
