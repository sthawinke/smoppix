#' Estimate the PI for the distance to a fixed region of interest
#'
#' @param pSub The subset point pattern containing only a single gene
#' @param ecdfAll the cumulative distribution function under the null
#' @inheritParams estPimsSingle
#' @param pi The type of PI to calculate
#'
#' @return The estimated probabilistic index
#' @importFrom spatstat.geom nncross edges nndist
#' @importFrom stats dist
#' @importFrom Rfast rowMins colMins
#' @details Analysis of the distance to the border was introduced by \insertCite{Joyner2013}{spatrans} in the form of the B-function.
#' The independent evaluations of the B-functions per cell are here returned as realizations of the probabilistic index.
#' @seealso \link{addCell}, \link{estPims}
#' @references
#' \insertAllCited{}
calcWindowDistPI = function(pSub, owins, centroids, ecdfAll, pi, null, ecdfs){
    splitPPP = split.ppp(pSub, f = "cell")
    obsDistEdge = lapply(names(splitPPP), function(x){
        Dist = switch(pi,
            "midpoint" = crossdist(splitPPP[[x]], centroids[[x]]),
            "edge" = nncross(splitPPP[[x]], edges(owins[[x]]), what = "dist"),
            "nnCell" = nndist(splitPPP[[x]]),
            "allDistCell" = dist(coords(splitPPP[[x]])))
        if(pi %in% c("edge", "midpoint")){
            ecdfAll[[x]][[pi]](Dist) #Do not average here, independent observations
        } else {
            if((NP <- npoints(splitPPP[[x]])) == 1)
                return(NA)
            if(pi == "allDistCell"){
                switch(null,
                   "CSR" = mean(ecdfAll[[x]]$allDistCell(Dist)),
                   "background" = {
                       Dist = as.matrix(Dist)
                       mean(vapply(seq_len(ncol(Dist)), FUN.VALUE = double(ncol(Dist)-1),
                                   function(i){
                                       ecdfs[[x]]$allDistCell[[i]](Dist[i, -i])
                                       }))
                       }
                       )
            } else if(pi == "nnCell"){
            approxRanks = switch(null,
                "CSR" = getApproxRanks(ecdfAll[[x]]$allDistCell, Dist),
                "background" = vapply(seq_along(Dist), FUN.VALUE = double(1),
                        function(i){getApproxRanks(ecdfs[[x]]$allDistCell[[i]], Dist[i])})
            )
            mean(vapply(seq_along(approxRanks), FUN.VALUE = double(1), function(ar){
                pnhyper(approxRanks[ar], r = 1, m = NP - 1,
                    n = switch(null, "CSR" = getN(ecdfAll[[x]]$allDistCell),
                               "background" = getN(ecdfs[[ar]]) - NP + 1))
            }))
        }
        }
    })
    names(obsDistEdge) = names(splitPPP)
    obsDistEdge
}
calcWindowPairPI = function(pSub1, pSub2, feat1, feat2, cd, ecdfAll, pi, null, ecdfs){
    splitPPP1 = split.ppp(pSub1, f = "cell")
    splitPPP2 = split.ppp(pSub2, f = "cell")
    out = vapply(nam <- intersect(names(splitPPP1), names(splitPPP2)),
                 FUN.VALUE = double(1), function(x){
        cd = crossdist(splitPPP1[[x]], splitPPP2[[x]])
        id1 = marks(splitPPP1[[x]], drop = FALSE)$gene == feat1
        id2 = marks(splitPPP2[[x]], drop = FALSE)$gene == feat2
        if(pi == "allDistPairCell"){
            switch(null,
                   "CSR" =  mean(ecdfAll[[x]]$allDistCell(cd)),
                   "background" = {
                       pi1 = vapply(seq_len(NP1 <- npoints(splitPPP1[[x]])),
                                    FUN.VALUE = double(ncol(cd)),
                                   function(i){ecdfs[[x]]$allDistCell[[id1[i]]](cd[i, ])})
                       pi2 = vapply(seq_len(npoints(splitPPP2[[x]])),
                                    FUN.VALUE = double(nrow(cd)),
                                    function(i){ecdfs[[x]]$allDistCell[[id2[i]]](cd[,i])})
                       mean(c(pi1, pi2))
                   }
            )
        } else if(pi == "nnPairCell"){
            nnDist = c(rowMins(cd, value = TRUE), colMins(cd, value = TRUE))
            np1 = npoints(splitPPP1[[x]]); np2 = npoints(splitPPP2[[x]])
            seq1 = seq_len(np1)
            if(null == "CSR"){
                approxRanks = getApproxRanks(ecdfAll[[x]]$allDistCell, nnDist)
                n = getN(ecdfAll[[x]]$allDistCell)
                pi1 = pnhyper(approxRanks[seq1], r = 1, m = np1, n = n - np1)
                pi2 = pnhyper(approxRanks[-seq1], r = 1, m = np2, n = n - np2)
            } else if(null == "background"){
                approxRanks = vapply(seq_along(nnDist), FUN.VALUE = double(1), function(i){
                    getApproxRanks(ecdfs[[x]]$allDistCell[[i]], nnDist[i])
                })
                pi1 = vapply(seq_along(approxRanks)[seq1], FUN.VALUE = double(1), function(ar){
                    pnhyper(approxRanks[ar], r = 1, m = np1 - 1,
                            n = getN(ecdfs[[x]]$allDistCell[[ar]]) - np1 + 1)
                })
                pi2 = vapply(seq_along(approxRanks)[-seq1], FUN.VALUE = double(1), function(ar){
                    pnhyper(approxRanks[ar], r = 1, m = np2 - 1,
                            n = getN(ecdfs[[x]]$allDistCell[[ar]]) - np2 + 1)
                })
            }
            mean(c(pi1, pi2))
        }
    })
    names(out) = nam
    return(out)
}
#' Approximate the ranks given and ecdf
#'
#' @param ecdf An empirical density function
#' @param obs A vector of observations
#'
#' @return The approximate ranks
getApproxRanks = function(ecdf, obs){
    ranks = round(ecdf(obs)*environment(ecdf)$nobs)
    ranks[ranks==0] = 1
    ranks
}

