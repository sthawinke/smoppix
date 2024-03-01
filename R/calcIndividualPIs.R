#' Calculate individual PI entries
#'
#' @param pSubLeft The subsampled overall point pattern returned by subSampleP
#' @param ecdfAll,ecdfsCell Empirical cumulative distribution functions of all
#' events and of cells within the cell, under the null
#' @inheritParams estPisSingle
#'
#' @return A list containing PI entries per feature
#' @details For the single-feature nearest neighbour distances, the average
#' is already calculated
#' @importFrom spatstat.geom nndist npoints nncross.ppp
calcIndividualPIs <- function(p, tabObs, pis, pSubLeft, owins, centroids, null,
                              features, ecdfAll, ecdfsCell, loopFun, minDiff, minObsNN) {
    NPall <- switch(null,
        "CSR" = max(tabObs) * 4,
        # NPall is arbitrary for CSR, make sure it is large enough,
        # while limiting computation time in evaluating the negative
        # hypergreometric
        "background" = npoints(p)
    )
    pSplit <- split.ppp(p, f = factor(marks(p, drop = FALSE)$gene))
    # Divide the work over the available workers
    piList <- loadBalanceBplapply(loopFun = loopFun, iterator = features, func = function(feat) {
        pSub <- pSplit[[feat]]
        if((NP <- npoints(pSub)) >= minObsNN){
            if (null == "background" && ((NPall - NP) < minDiff)) {
                distMat <- NULL
                calcNNsingle <- FALSE
                # If insufficient other events, do not calculate nn PIs for background
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
                    # Cross-distances
                })
            }

            if (isMat <- is.matrix(distMat)) {
                approxRanksTmp <- switch(null,
                    "background" = findRanksDist(
                        getCoordsMat(pSub), getCoordsMat(pSubLeft$Pout),
                        distMat^2
                    ),
                    "CSR" = matrix(ecdfAll(distMat), nrow = nrow(distMat))
                )
                approxRanks <- round((approxRanksTmp/(switch(null,
                        "background" = (npoints(pSubLeft$Pout) -
                                          which(marks(p, drop = FALSE)$gene == feat)
                                      %in% pSubLeft$id),
                        "CSR" = 1
                ))) * NPall)
                # Correct for self distances by subtracting one. The C++ function only counts distances larger so
                # self distances will be zero in the numerator
                #The observed nearest neighbour distances remain part of the permutation,
                #See Smyth "permutation p-values can never be zero"
                #Need binomial here, and stop halfway the quantile in the negative hypergeometric?
                # Get the order right to prevent integer overflow: first divide, then multiply
                approxRanks[approxRanks == 0] <- 1
                colnames(approxRanks) <- colnames(distMat)
            }
            # Names may get lost in C++ function And then rearrange to get to the PIs
            nnPI <- if (calcNNsingle && isMat) {
                mean(calcNNPI(approxRanks[, "self"], NPall - (NP - 1), m = NP - 1, r = 1))
            } else {
                NA
            }
            nnPIpair <- if ("nnPair" %in% pis && isMat) {
                vapply(setdiff(colnames(approxRanks), "self"), FUN.VALUE = double(NP), function(g) {
                    NP <- tabObs[g]
                    calcNNPI(approxRanks[, g], n = NPall - NP, m = NP, r = 1)
                })
            }
        } else {
            nnPI = nnPIpair = NA
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
