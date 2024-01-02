#' Estimate the PI for the nearest neighbour distances
#'
#' @inheritParams estPimsSingle
#' @param pSub The subset point pattern with only one feature
#' @param ecdfAll The empirical cumulative distribution function of
#' all distances under the null
#' @param n the number of subsampled distances
#'
#' @importFrom spatstat.geom nndist crossdist npoints area
#' @importFrom spatstat.random rpoispp
#' @importFrom stats ecdf
#' @importFrom Rfast rowMins colMins
#' @importFrom extraDistr pnhyper
#' @return The estimated probabilistic index
calcNNPI <- function(null, id, p, ecdfAll) {
    NP <- length(id)
    if(NP <= 1 ){
        return(NA)
    } else {
        obsDistNN <- nndist(p[id,])
    }
    if (null == "background") {
        nFac = n/ncol(cd)
        approxRanks <- round(nFac*(rowSumsWrapper(cd, obsDistNN) + 0.5))
        mean(pnhyper(approxRanks, n = n - (NP - 1), m = NP - 1, r = 1))
        # Exclude event itself, so NP - 1 m = N-1: White balls, number of other
        #events of the same gene n = n - (NP-1): black balls, number of events
        #of other genes in background r=1: Nearest neighbour so first
        # occurrence
    } else if (null == "CSR") {
        # Weigh by negative hypergeometric distribution to
        #bypass Monte-Carlo simulations
        approxRanks <- getApproxRanks(ecdfAll, obsDistNN, n)
        mean(pnhyper(approxRanks, n = n - (NP - 1), m = NP - 1, r = 1))
    }
}
calcNNPIpair <- function(obsDistNN, id1, id2, null, cd, n, ecdfAll) {
    obsDistRank <- if (null == "background") {
        nFac = n/ncol(cd)
        approxRanks <- round(nFac*(rowSumsWrapper(cd, obsDistNN) + 0.5))
        approxRanks
    } else {
        getApproxRanks(ecdfAll, obsDistNN, n)
    }
    seq1 <- seq_along(id1)
    pis <- mean(c(pnhyper(obsDistRank[seq1], n = n - length(id2),
                          m = length(id2), r = 1),
                  pnhyper(obsDistRank[-seq1], n = n - length(id1),
                          m = length(id1), r = 1)))
    # Using the negative hypergeometric precludes Monte-Carlo
    return(pis)
}
