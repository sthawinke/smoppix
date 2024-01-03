#' Calculate individual PI entries
#'
#' @param pSubLeft The subsampled overall point pattern
#' @param ecdfAll,ecdfsCell Empirical cumulative distribution functions of all
#' events and of cells within the cell, under the null
#' @inheritParams estPimsSingle
#'
#' @return A list over the features with all PI entries
#' @details For the single-feature nearest neighbour distances, the average
#' is already calculated
#' @importFrom spatstat.random runifpoint
#' @importFrom utils combn
#' @importFrom spatstat.geom nncross crossdist coords npoints
#' @importFrom BiocParallel bpparam bplapply
calcIndividualPIs = function(p, tabObs, pis, pSubLeft, owins, centroids, null,
                             features, ecdfAll, ecdfsCell){
    NPall = npoints(p)
    pSplit = split.ppp(p, f = factor(marks(p, drop = FALSE)$gene))
    splitFac <- rep(seq_len(bpparam()$workers), length.out = length(nams <- names(tabObs[features])))
    #Divide the work over the available workers
    piList <- unsplit(f = splitFac, bplapply(split(nams, f = splitFac), function(ss){
        featPIs = lapply(ss, function(feat) {
            pSub <- pSplit[[feat]]
            NP <- npoints(pSub)
            distMat = cbind(
                "self" = if(calcNNsingle <- (NP> 1 && ("nn" %in% pis))){
                    nndist(pSub)}, #Self distances
                if("nnPair" %in% pis){
                    id = !(names(pSplit) %in% feat)
                    matrix(unlist(lapply(pSplit[id],
                                         function(y){nncross(pSub, y, what = "dist")
                                         })), nrow = NP, dimnames = list(NULL, names(pSplit)[id]))
                    #Cross-distances
                })
            if(is.matrix(distMat)){
                approxRanksTmp = switch(null,
                    "background" = findRanksDist(getCoordsMat(pSub),
                                                 getCoordsMat(pSubLeft$Pout), distMat^2),
                    "CSR" = matrix(ecdfAll(distMat), nrow = nrow(distMat))
                )
                approxRanks = round(NPall*approxRanksTmp/
                    (switch(null,
                            "background" = (npoints(pSubLeft$Pout)-
                                                                                     which(marks(p, drop = FALSE)$gene == feat) %in% pSubLeft$id),
                                                                 #Correct for self distances by subtracting one.
                                                                 #The C++ function only counts distances larger so self distances
                                                                 #will be ignored in the numerator
                                                                 "CSR" = 1)))
                approxRanks[approxRanks==0] = 1
                colnames(approxRanks) = colnames(distMat)
            }
            #Names may get lost in C++ function
            #And then rearrange to get to the PIs
            nnPI = if(calcNNsingle){
                mean(calcNNPI(approxRanks[, "self"], NPall - (NP - 1),
                              m = NP - 1, r = 1))
            } else NA
            nnPIpair = if("nnPair" %in% pis){
                apply(approxRanks[, colnames(approxRanks) != "self", drop = FALSE], 2,
                      calcNNPI, n = NPall - NP, m = NP, r = 1)
            }
            ## Window related distances
            edgeDistPI <- if (any(pis == "edge")) {
                calcWindowDistPI(pSub, owins = owins, ecdfAll = ecdfsCell,
                                 pi = "edge")
            }
            midPointDistPI <- if (any(pis == "midpoint")) {
                calcWindowDistPI(pSub, centroids = centroids,
                                 ecdfAll = ecdfsCell, pi = "midpoint")
            }
            list("windowDists" = list("edge" = edgeDistPI,
                                      "midpoint" = midPointDistPI),
                 "pointDists" = list("nn" = nnPI, "nnPair" = nnPIpair))
        })
        names(featPIs) = ss
        # TO DO: nnCell recycling code
        return(featPIs)
    }))
    names(piList) = nams
    return(piList)
}