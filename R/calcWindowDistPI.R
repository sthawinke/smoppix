#' Estimate the PI for the distance to a fixed region of interest
#'
#' @param pSub The subset point pattern containing only a single gene
#' @param ecdfAll the cumulative distribution function under the null
#' @param owins the list of windows describing the cells
#' @param pi a character string, the type of PI desired
#'
#' @return The estimated probabilistic index
#' @importFrom spatstat.geom nncross centroid.owin edges nndist
#' @importFrom stats dist
#' @importFrom Rfast rowMins colMins
#' @details Analysis of the distance to the border was introduced by \insertCite{Joyner2013}{spatrans} in the form of the B-function.
#' The independent evaluations of the B-functions per cell are here returned as realizations of the probabilistic index.
#' @seealso \link{addCell}, \link{estPims}
#' @references
#' \insertAllCited{}
calcWindowDistPI = function(pSub, owins, ecdfAll, pi){
    splitPPP = split.ppp(pSub, f = "cell")
    obsDistEdge = lapply(names(splitPPP), function(x){
        Dist = switch(pi,
            "edge" = crossdist(splitPPP[[x]], centroid.owin(owins[[x]], as.ppp = TRUE)),
            "midpoint" = nncross(splitPPP[[x]], edges(owins[[x]]), what = "dist"),
            "nnCell" = nndist(splitPPP[[x]]),
            "allDistCell" = dist(coords(splitPPP[[x]])))
        if(pi %in% c("edge", "midpoint")){
            ecdfAll[[x]][[pi]](Dist) #Do not average here, independent observations
        } else if(pi == "allDistCell"){
            mean(ecdfAll[[x]]$allDistCell(Dist))
        } else if(pi == "nnCell"){
            approxRanks = getApproxRanks(ecdfAll[[x]]$allDistCell, Dist)
            getPoissonPi(npoints(splitPPP[[x]]),
                         nAll = environment(ecdfAll[[x]]$allDistCell)$nobs, approxRanks)
        }
    })
    names(obsDistEdge) = names(splitPPP)
    obsDistEdge
}
calcWindowPairPI = function(pSub1, pSub2, cd, ecdfAll, pi){
    splitPPP1 = split.ppp(pSub1, f = "cell")
    splitPPP2 = split.ppp(pSub2, f = "cell")
    lapply(intersect(names(splitPPP1), names(splitPPP2)), function(x){
        cd = crossdist(splitPPP1[[x]], splitPPP2[[x]])
        if(pi == "allDistPairCell"){
            mean(ecdfAll[[x]]$allDistCell(cd))
        } else if(pi == "nnPairCell"){
            nnDist = c(rowMins(cd), colMins(cd))
            approxRanks = getApproxRanks(ecdfAll[[x]]$allDistCell, nnDist)
            nAll = environment(ecdfAll[[x]]$allDistCell)$nobs
            seq1 = seq_len(nrow(cd))
            pi1 = getPoissonPi(np1, nAll, approxRanks[seq1])
            pi2 = getPoissonPi(np2, nAll, approxRanks[-seq1])
            (pi1*np1+pi2*np2)/(np1+np2)
        }
    })
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

