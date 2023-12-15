#' Estimate the PI for the nearest neighbour distances
#'
#' @inheritParams estPimsSingle
#' @inheritParams calcAllDistPI
#' @param n the number of subsampled distances
#'
#' @importFrom spatstat.geom nndist crossdist npoints area
#' @importFrom spatstat.random rpoispp
#' @importFrom stats ecdf
#' @importFrom Rfast rowMins colMins
#' @importFrom extraDistr pnhyper
#' @return The estimated probabilistic index
calcNNPI <- function(pSub, p, null, ecdfs, n, ecdfAll) {
    obsDistNN <- nndist(pSub)
    NP <- npoints(pSub)
    if (null == "background") {
        approxRanks <- vapply(seq_along(obsDistNN),
            FUN.VALUE = double(1),
            function(i) {
                round(ecdfs[[i]](obsDistNN[i]) * n)
                # Approximate rank: quantile in overall distribution times
                # number of observations.
            }
        )
        approxRanks[approxRanks == 0] <- 1
        mean(pnhyper(approxRanks, n = n - (NP - 1), m = NP - 1, r = 1))
        # Exclude event itself, so NP - 1
        # m = N-1: White balls, number of other events of the same gene
        # n = n - (NP-1): black balls, number of events of other genes in background
        # r=1: Nearest neighbour so first occurrence
    } else if (null == "CSR") {
        # Weigh by Poisson and negative hypergeometric distribution to bypass
        # Monte-Carlo simulations
        approxRanks <- getApproxRanks(ecdfAll, obsDistNN)
        mean(pnhyper(approxRanks, n = getN(ecdfAll) - (NP - 1), m = NP - 1, r = 1))
    }
}
calcNNPIpair <- function(obsDistNN, id1, id2, null, p, ecdfs, n, ecdfAll) {
    obsDistRank <- if (null == "background") {
        vapply(seq_along(obsDistNN), FUN.VALUE = double(1), function(i) {
            round(ecdfs[[i]](obsDistNN[i]) * n)
        })
    } else {
        round(ecdfAll(obsDistNN) * (n <- getN(ecdfAll)))
    }
    obsDistRank[obsDistRank == 0] <- 1
    seq1 <- seq_along(id1)
    pis <- mean(c(
        pnhyper(obsDistRank[seq1], n = n - length(id2), m = length(id2), r = 1),
        pnhyper(obsDistRank[-seq1], n = n - length(id1), m = length(id1), r = 1)
    ))
    # Using the negative hypergeometric precludes Monte-Carlo
    return(pis)
}
