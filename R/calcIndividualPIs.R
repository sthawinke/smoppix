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
#' @seealso \link{estPis}, \link{calcNNPI}
calcIndividualPIs <- function(p, tabObs, pis, pSubLeft, owins, centroids, null,
                              features, ecdfAll, ecdfsCell, loopFun, minDiff, minObsNN) {
    NPall <- if(bg <- (null=="background")){
        npoints(p)
    } else if(null=="CSR"){
        max(tabObs) * 4
        # NPall is arbitrary for CSR, make sure it is large enough,
        # while limiting computation time in evaluating the negative
        # hypergeometric
    }
    pSplit <- split.ppp(p, f = factor(marks(p, drop = FALSE)$gene))
    # Divide the work over the available workers
    piList <- loadBalanceBplapply(loopFun = loopFun, iterator = features, func = function(feat) {
        pSub <- pSplit[[feat]]
        if((NP <- npoints(pSub)) >= minObsNN){
            if (bg && ((NPall - NP) < minDiff)) {
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
                })
            }
            if (isMat <- is.matrix(distMat)) {
                featId = which(marks(p, drop = FALSE)$gene == feat)
                # cd = crossdistWrapper(pSub, pSubLeft$Pout)
                # approxRanksNew = rowSums(distMat[, "self"] > cd)
                # tiesMatNew = rowSums(distMat[, "self"] == cd)
                doubleMatrixRanks <- switch(null,
                    "background" = findRanksDist(
                        getCoordsMat(pSub), getCoordsMat(pSubLeft$Pout),
                        distMat
                    ),
                    "CSR" = matrix(ecdfAll(distMat), nrow = nrow(distMat))
                )
                approxRanksTmp = doubleMatrixRanks[seq_len(nrow(distMat)),,drop = FALSE]
                tiesMat = doubleMatrixRanks[-seq_len(nrow(distMat)),,drop = FALSE]
                #Leave out found nearest neighbour distance in calculation,
                #and add 1 to ties
                if(bg) {
                    selfPoint = (featId %in% pSubLeft$id)
                    selfPointRanks = selfPoint & distMat > 0
                    selfPointTies = selfPoint & distMat == 0
                    approxRanks <- round(((approxRanksTmp-selfPointRanks)/(npoints(pSubLeft$Pout) -selfPoint)
                    ) * (NPall-selfPoint))
                    tiesMat <- round(((tiesMat-selfPointTies)/(npoints(pSubLeft$Pout) - selfPoint)
                    ) * (NPall-selfPoint))
                    tiesMat[tiesMat==0] = 1
                    #Cfr permutation p-values can never be zero (Phipson 2010)
                    #The observed distance is at least a tie.
                    #Adding 1 to the denominator too will not make a big difference
                } else {
                    approxRanks <- round(approxRanksTmp*NPall)
                }

                # Correct for self distances by subtracting one everywhere.
                #The observed nearest neighbour distances (not the same as self distances) remain part of the permutation,
                # Get the order right to prevent integer overflow: first divide, then multiply
                colnames(approxRanks) <- colnames(tiesMat) <- colnames(distMat)
            }
            # Names may get lost in C++ function. Then rearrange to get to the PIs
            nnPI <- if (calcNNsingle && isMat) {
                mean(calcNNPI(approxRanks[, "self"], NPall - (NP - 1) - bg,
                              m = NP - 1, r = 1, ties = if(bg) tiesMat[, "self"]))
            } else {
                NA
            }
            #Minus bg subtracts self distances
            nnPIpair <- if ("nnPair" %in% pis && isMat) {
                vapply(setdiff(colnames(approxRanks), "self"), FUN.VALUE = double(NP), function(g) {
                    NP <- tabObs[g] #Number of other gene in the pair
                    calcNNPI(approxRanks[, g], n = NPall - NP - bg, m = NP,
                             r = 1, ties = if(bg) tiesMat[, g])
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
