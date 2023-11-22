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
#' @return The estimated probabilistic index
calcNNPI = function(pSub, p, null, nSims, ecdfs, n){
    obsDistNN = nndist(pSub)
    if(null == "background"){
        NP = npoints(pSub)
        obsDistRank = vapply(seq_along(obsDistNN), FUN.VALUE = double(1), function(i){
                round(ecdfs[[i]](obsDistNN[i])*environment(ecdfs[[i]])$nobs)
            }) + 1
        mean(pnhyper(obsDistRank, n = n- NP + 1, m = NP - 1, r = 1))#Exclude event itself
    } else if(null == "CSR"){
        lambda = npoints(pSub)/area(p$window)
        simDistsNN = lapply(integer(nSims), function(i){
            nndist(rpoispp(lambda, win = pSub$window))
        }) #Density dependent so still sims needed
        mean(ecdf(unlist(simDistsNN))(obsDistNN))
    }
}
calcNNPIpair = function(cd, id1, id2, null, p, ecdfs, nSims){
    npp = npoints(p)
    obsDistNN = c(rowMins(cd, value = TRUE), colMins(cd, value = TRUE))
    if(null == "background"){
        obsDistRank = vapply(c(id1, id2), FUN.VALUE = double(1), function(i){
            round(ecdfs[[i]](obsDistNN[i])*environment(ecdfs[[i]])$nobs)
        }) + 1
        seq1 = seq_along(id1)
        PP = ncol(obsDistRank)
        mean(c(pnhyper(obsDistRank[seq1], n = PP-length(id2), m = length(id2), r = 1),
                             pnhyper(obsDistRank[-seq1], n = PP-length(id1), m = length(id1), r = 1)))
        #Using the negative hypergeometric precludes Monte-Carlo
    } else if(null == "CSR"){
        pSim = runifpoint(npp, win = p$window)
        np1 = length(id1);np2 = length(id2)
        simDistsNN = vapply(FUN.VALUE = obsDistNN, integer(nSims), function(i){
            id <- sample(npp, np1+np2)
            p1 = pSim[id[seq_len(np1)],];p2 = pSim[id[seq_len(np2)],]
            distMatSub = crossdist(p1, p2)
            c(rowMins(distMatSub, value = TRUE), colMins(distMatSub, value = TRUE))
        })
        mean(simDistsNN < obsDistNN)
    }
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
