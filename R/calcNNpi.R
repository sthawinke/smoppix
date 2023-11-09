#' Estimate the PI for the nearest neighbour distances
#'
#' @param pSub The subset point pattern containing only a single gene
#' @inheritParams estPimsSingle
#'
#' @importFrom matrixStats rowMins
#' @importFrom spatstat.geom nndist crossdist npoints area
#' @importFrom spatstat.random rpoispp
#' @importFrom stats ecdf
#' @return The estimated probabilistic index
calcNNPI = function(pSub, p, null, nSims){
    obsDistNN = nndist(pSub)
    if(null == "background"){
        simDistsNN = vapply(integer(nSims), FUN.VALUE = double(npoints(pSub)), function(i){
            distMat = crossdist(pSub, subSampleP(p, npoints(pSub)-1))
            # #Keep observed points fixed
            matrixStats::rowMins(distMat)
            #Keep this for bivariate analysis: Average over all points, single ecdf per point
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
calcNNPIpair = function(pSub1, pSub2, p, null, nSims, distMat){
    obsDistNN1 = nncross(pSub1, pSub2, what = "dist");obsDistNN2 = nncross(pSub2, pSub1, what = "dist")
    np1 = npoints(pSub1);np2 = npoints(pSub2); npTot <- (np1 + np2)
    if(null == "background"){
        simDistsNN1 = vapply(integer(nSims), FUN.VALUE = double(np1), function(i){
            nncross(pSub1, subSampleP(p, np2), what = "dist")
        })
        simDistsNN2 = vapply(integer(nSims), FUN.VALUE = double(np2), function(i){
            nncross(pSub2, subSampleP(p, np1), what = "dist")
        })
        isMat1 = np1 > 1;isMat2 = np2 > 1
        piEsts1 = vapply(FUN.VALUE = double(1), seq_len(np1), function(i){
            ecdf(if(isMat1) simDistsNN1[i,] else simDistsNN1)(obsDistNN1[i])
        })
        piEsts2 = vapply(FUN.VALUE = double(1), seq_len(np2), function(i){
            ecdf(if(isMat2) simDistsNN2[i,] else simDistsNN2)(obsDistNN2[i ])
        })
        mean(c(piEsts1, piEsts2))
    } else if(null == "CSR"){
        simDistsNN = lapply(integer(nSims), function(i){
            pSim = runifpoint(npTot, win = p$window)
            p1 = pSim[id <- sample(npTot, np1),];p2 = pSim[-id,]
            c(nncross(what = "dist", p1, p2), nncross(what = "dist", p2, p1))
        })
        mean(ecdf(unlist(simDistsNN))(c(obsDistNN1, obsDistNN2)))
    }
}
