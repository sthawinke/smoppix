#' Calculate individual PI entries of a single point pattern
#'
#' @param pSubLeft The subsampled overall point pattern returned by subSampleP
#' @param ecdfAll,ecdfsCell Empirical cumulative distribution functions of all
#' events and of cells within the cell, under the null
#' @inheritParams estPisSingle
#'
#' @return A list containing PI entries per feature
#' @details For the single-feature nearest neighbour distances, the PI is average
#' over the point pattern
#' @importFrom spatstat.geom nndist npoints nncross.ppp
#' @seealso \link{estPis}, \link{calcNNPI}
calcIndividualPIs <- function(p, tabObs, pis, pSubLeft, owins, centroids, null,
                              features, ecdfAll, ecdfsCell, loopFun, minDiff, minObsNN) {
    NPall <- if (bg <- (null == "background")) {
        npoints(p)
    } else if (null == "CSR") {
        max(tabObs) * 4
        # NPall is arbitrary for CSR, make sure it is large enough,
        # while limiting computation time in evaluating the negative
        # hypergeometric
    }
    pSplit <- split.ppp(p, f = factor(marks(p, drop = FALSE)$gene))
    # Divide the work over the available workers
    piList <- loadBalanceBplapply(loopFun = loopFun, iterator = features, func = function(feat) {
        pSub <- pSplit[[feat]] #ppp-object of a single feature
        if ((NP <- npoints(pSub)) >= minObsNN) {
            if (bg && ((NPall - NP) < minDiff)) {
                distMat <- NULL
                calcNNsingle <- FALSE
                # If insufficient number of other events, do not calculate nn PIs for background
            } else {
                distMat <- cbind(self = if (calcNNsingle <- (NP > 1 && ("nn" %in% pis))) {
                    nndist(pSub)
                }, if ("nnPair" %in% pis) {
                    id <- !(names(pSplit) == feat)
                    if (any(id)) {
                        matrix(unlist(lapply(pSplit[id], function(y) {
                            nncross.ppp(pSub, y,
                                what = "dist", sortby = "x",
                                is.sorted.X = TRUE, is.sorted.Y = TRUE
                            )
                            # Point patterns have been pre-sorted in hyperframe function
                        })), nrow = NP, dimnames = list(NULL, names(pSplit)[id]))
                    }
                }) #Observed distances
            }
            if (isMat <- is.matrix(distMat)) {
                doubleMatrixRanks <- switch(null,
                    "background" = findRanksDist(
                        getCoordsMat(pSub), getCoordsMat(pSubLeft$Pout),
                        distMat
                    ),
                    "CSR" = matrix(ecdfAll(distMat), nrow = nrow(distMat))
                )
                approxRanks <- doubleMatrixRanks[seq_len(nrow(distMat)), , drop = FALSE]
                tiesMat <- doubleMatrixRanks[-seq_len(nrow(distMat)), , drop = FALSE]
                if (bg) {
                   featId <- match(feat, marks(p, drop = FALSE)$gene)
                   #Indices of feature in subsampled ppp
                    selfPoint <- (featId %in% pSubLeft$id)
                    # Correct for self distances, when point itself is part of the permutation, 
                    # leading to distance 0, by subtracting one everywhere.
                    # The observed nearest neighbour distances (not the same as self distances) remain part of the permutation
                    selfPointRanks <- selfPoint & distMat > 0
                    selfPointTies <- selfPoint & distMat == 0
                    approxRanks <- round(((approxRanks - selfPointRanks) / (npoints(pSubLeft$Pout) - selfPoint)
                    ) * (NPall - selfPoint))
                    # Get the order right to prevent integer overflow: first divide, then multiply
                    tiesMatCor <- tiesMat - selfPointTies
                    tiesMatCorId <- tiesMatCor == 0
                    # Cfr permutation p-values can never be zero (Phipson 2010)
                    # The observed distance is at least a tie.
                    tiesMat <- round(((tiesMatCor+tiesMatCorId) / (npoints(pSubLeft$Pout) - selfPoint + tiesMatCorId)
                    ) * (NPall - selfPoint))
                  } else {
                    approxRanks <- round(approxRanks * NPall)
                }
                colnames(approxRanks) <- colnames(tiesMat) <- colnames(distMat)
                # Names get lost in C++ function.
            }
            nnPI <- if (calcNNsingle && isMat) {
                mean(calcNNPI(approxRanks[, "self"], NPall - (NP - 1) - bg,
                    m = NP - 1, r = 1, ties = if (bg) tiesMat[, "self"]
                ))
            } else {
                NA
            }
            # Minus bg subtracts self distances
            nnPIpair <- if ("nnPair" %in% pis && isMat) {
                vapply(setdiff(colnames(approxRanks), "self"), FUN.VALUE = double(NP), function(g) {
                    NP <- tabObs[g] # Number of other gene in the pair
                    calcNNPI(approxRanks[, g],
                        n = NPall - NP - bg, m = NP,
                        r = 1, ties = if (bg) tiesMat[, g]
                    )
                })
            }
        } else {
            nnPI <- nnPIpair <- NA
        }
        ## Window related distances
        edgeDistPI <- if (any(pis == "edge")) {
            calcWindowDistPI(pSub, owins = owins, ecdfAll = ecdfsCell, pi = "edge")
        }
        midPointDistPI <- if (any(pis == "centroid")) {
            calcWindowDistPI(pSub, centroids = centroids, ecdfAll = ecdfsCell, pi = "centroid")
        }
        list(windowDists = list(edge = edgeDistPI, centroid = midPointDistPI), pointDists = list(nn = nnPI, nnPair = nnPIpair))
    })
    names(piList) <- features
    return(piList)
}
