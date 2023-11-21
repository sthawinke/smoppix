#' Estimate the PI for the nearest neighbour distances
#'
#' @param pSub The subset point pattern containing only a single gene
#' @inheritParams estPimsSingle
#'
#' @importFrom spatstat.geom nndist crossdist npoints area
#' @importFrom spatstat.random rpoispp
#' @importFrom stats ecdf
#' @importFrom Rfast rowMins colMins
#' @return The estimated probabilistic index
calcNNPI = function(pSub, p, null, nSims){
    obsDistNN = nndist(pSub)
    if(null == "background"){
        simDistsNN = vapply(integer(nSims), FUN.VALUE = double(npoints(pSub)), function(i){
            nncrossFast(pSub, subSampleP(p, npoints(pSub)-1))
            # #Keep observed points fixed
        })
        mean(simDistsNN < obsDistNN)
    } else if(null == "CSR"){
        lambda = npoints(pSub)/area(p$window)
        simDistsNN = lapply(integer(nSims), function(i){
            nndist(rpoispp(lambda, win = pSub$window))
        })
        mean(ecdf(unlist(simDistsNN))(obsDistNN))
    }
}
calcNNPIpair = function(pSub1, pSub2, p, pJoin, null, nSims, ecdfsNN, feat1, feat2, tabObs){
    distMatSub = crossdist(pSub1, pSub2)
    obsDistNN1 = rowMins(distMatSub, value = TRUE)
    obsDistNN2 = colMins(distMatSub, value = TRUE)
    # np1 = npoints(pSub1);np2 = npoints(pSub2); npTot <- (np1 + np2);npp = npoints(p)
    # Replace = npTot > npp
    if(null == "background"){
        # tmpDist = crossdist(pJoin, p)
        # #This breaks the dependence, but we are averaging anyway.
        # #The dependence between distances merely inflates the variability of the PI, it does not affect its null distribution
        # sam = if(Replace) {
        #     rnbinom(nSims, size = 1, prob = 1 - (1-1/npTot)^npp)+1
        # } else {
        #     rnhyper(nSims, m = npTot, n = npp-npTot, r = 1)}
        # simDistsNN = rowSort(tmpDist)[, sam, drop = FALSE]
        # # simDistsNN = vapply(integer(nSims), FUN.VALUE = double(npTot), function(i){
        # #     rowMins(tmpDist[, sample(npp, npTot, replace = Replace), drop = FALSE])
        # # })
        mean(c(mapply(ecdfsNN[[feat1]][[as.character(tabObs[feat2])]], obsDistNN1, FUN = function(e, x) e(x)),
               mapply(ecdfsNN[[feat2]][[as.character(tabObs[feat1])]], obsDistNN2, FUN = function(e, x) e(x))))
    } else if(null == "CSR"){
        pSim = runifpoint(npp <- npoints(p), win = p$window)
        simDistsNN = lapply(integer(nSims), function(i){
            id <- sample(npp, npTot, replace = Replace)
            p1 = pSim[id[seq_len(np1)],];p2 = pSim[-id[seq_len(np1)],]
            distMatSub = crossdist(p1, p2)
            c(rowMins(distMatSub, value = TRUE), colMins(distMatSub, value = TRUE))
        })
        mean(ecdf(unlist(simDistsNN))(obsDistNN))
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
