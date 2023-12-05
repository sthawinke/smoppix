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
calcNNPI = function(pSub, p, null, ecdfs, n, ecdfAll){
    obsDistNN = nndist(pSub)
    if(null == "background"){
        NP = npoints(pSub)
        approxRanks = vapply(seq_along(obsDistNN), FUN.VALUE = double(1), function(i){
                round(ecdfs[[i]](obsDistNN[i])*environment(ecdfs[[i]])$nobs)
            #Approximate rank: quantile in overall distribution times
            #number of observations.
            })
        approxRanks[approxRanks==0] = 1
        mean(pnhyper(approxRanks, n = n - (NP - 1), m = NP - 1, r = 1))
        #Exclude event itself, so NP - 1
        #m = N-1: White balls, number of other events of the same gene
        #n = n - (NP-1): black balls, number of events of other genes in background
        #r=1: Nearest neighbour so first occurrence
    } else if(null == "CSR"){
        #Weigh by Poisson and negative hypergeometric distribution to bypass Monte-Carlo simulations
        approxRanks = round(ecdfAll(obsDistNN)*(nAll <- environment(ecdfAll)$nobs))
        approxRanks[approxRanks==0] = 1
        getPoissonPi(NP = npoints(pSub), nAll, approxRanks)
    }
}
calcNNPIpair = function(cd, id1, id2, null, p, ecdfs, n, ecdfAll){
    npp = npoints(p)
    obsDistNN = c(rowMins(cd, value = TRUE), colMins(cd, value = TRUE))
    obsDistRank = if(null == "background"){
         vapply(seq_along(obsDistNN), FUN.VALUE = double(1), function(i){
            round(ecdfs[[i]](obsDistNN[i])*environment(ecdfs[[i]])$nobs)
        })
    } else {round(ecdfAll(obsDistsNN)*environment(ecdfAll)$nobs)}
        obsDistRank[obsDistRank==0] = 1
        seq1 = seq_along(id1)
    pis = if(null=="background"){
        mean(c(pnhyper(obsDistRank[seq1], n = n-length(id2), m = length(id2), r = 1),
               pnhyper(obsDistRank[-seq1], n = n-length(id1), m = length(id1), r = 1)))
        #Using the negative hypergeometric precludes Monte-Carlo
    } else if(null == "CSR"){
        np1 = length(id1);np2 = length(id2)
        # simDistsNN = vapply(FUN.VALUE = obsDistNN, integer(nSims), function(i){
        #     id <- sample(npp, np1+np2)
        #     p1 = pSim[id[seq_len(np1)],];p2 = pSim[id[seq_len(np2)],]
        #     distMatSub = crossdist(p1, p2)
        #     c(rowMins(distMatSub, value = TRUE), colMins(distMatSub, value = TRUE))
        # })
        # mean(simDistsNN < obsDistNN)
        pi1 = getPoissonPi(np1, nAll, obsDistRank[seq_len(np1)])
        pi2 = getPoissonPi(np2, nAll, obsDistRank[-seq_len(np1)])
        (pi1*np1+pi2*np2)/(np1+np2)
    }
    return(pis)
}
#' Fast version of nncross
#'
#' @param p1,p2 The point patterns
#'
#' @return A vector of nearest neighbour distances
#' @importFrom Rfast rowMins
#' @importFrom spatstat.geom crossdist
nncrossFast = function(p1, p2){
    rowMins(crossdist(p1,p2), value = TRUE)
}
