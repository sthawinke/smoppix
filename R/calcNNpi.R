#' Estimate the PI for the nearest neighbour distances
#'
#' @param pSub The subset point pattern containing only a single gene
#' @inheritParams estPimsSingle
#'
#' @importFrom matrixStats rowMins colMins
#' @importFrom spatstat.geom nndist crossdist npoints area
#' @importFrom spatstat.random rpoispp
#' @importFrom stats ecdf
#' @return The estimated probabilistic index
calcNNPI = function(pSub, p, null, nSims){
    obsDistNN = nndist(pSub)
    if(null == "background"){
        simDistsNN = vapply(integer(nSims), FUN.VALUE = double(npoints(pSub)), function(i){
            nncrossFast(pSub, subSampleP(p, npoints(pSub)-1))
            # #Keep observed points fixed
        })
        piEsts = vapply(seq_len(npoints(pSub)), FUN.VALUE = double(1), function(i){
            ecdf(simDistsNN[i,])(obsDistNN[i])
        })
        mean(piEsts)
    } else if(null == "CSR"){
        lambda = npoints(pSub)/area(p$window)
        simDistsNN = lapply(integer(nSims), function(i){
            nndist(rpoispp(lambda, win = pSub$window))
        })
        mean(ecdf(unlist(simDistsNN))(obsDistNN))
    }
}
calcNNPIpair = function(pSub1, pSub2, p, pJoin, null, nSims, distMat){
    distMatSub = crossdist(pSub1, pSub2)
    obsDistNN = c(rowMins(distMatSub), colMins(distMatSub))
    np1 = npoints(pSub1);np2 = npoints(pSub2); npTot <- (np1 + np2);npp = npoints(p)
    Replace = npTot > npp
    if(null == "background"){
        tmpDist = crossdist(pJoin, p)
        simDistsNN = vapply(integer(nSims), FUN.VALUE = double(npTot), function(i){
            rowMins(tmpDist[, sample(npp, npTot, replace = Replace), drop = FALSE])
        })
        piEsts = vapply(FUN.VALUE = double(1), seq_len(np1), function(i){
            ecdf(simDistsNN[i,])(obsDistNN[i])
        })
        mean(piEsts)
    } else if(null == "CSR"){
        pSim = runifpoint(npp <- npoints(p), win = p$window)
        simDistsNN = lapply(integer(nSims), function(i){
            id <- sample(npp, npTot, replace = Replace)
            p1 = pSim[id[seq_len(np1)],];p2 = pSim[-id[seq_len(np1)],]
            distMatSub = crossdist(p1, p2)
            c(rowMins(distMatSub), colMins(distMatSub))
        })
        mean(ecdf(unlist(simDistsNN))(obsDistNN))
    }
}
#' Fast version of nncross
#'
#' @param p1,p2 The point patterns
#'
#' @return A vector of nearest neighbour distances
#' @importFrom matrixStats rowMins
#' @importFrom spatstat.geom crossdist
nncrossFast = function(p1, p2){
    rowMins(crossdist(p1,p2))
}
