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
            #nncross(what = "dist", pSub, p[sample(npoints(p), npoints(p)),])
            distMat = crossdist(pSub, p[sample(npoints(p), npoints(pSub)-1),])
            distMat[distMat==0] = Inf
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
    obsDistNN = nncross(pSub1, pSub2, what = "dist")
    np1 = npoints(pSub1);np2 = npoints(pSub2); npTot <- (np1 + np2)
    if(null == "background"){
        simDistsNN = vapply(integer(nSims), FUN.VALUE = double(npTot), function(i){
            c(nncross(pSub1, p[sample(npoints(p), np2),], what = "dist"),
            nncross(pSub2, p[sample(npoints(p), np1),], what = "dist"))
        })
        piEsts = vapply(FUN.VALUE = double(1), seq_len(npTot), function(i){
            ecdf(simDistsNN[i,])(obsDistNN[i])
        })
        mean(piEsts)
    } else if(null == "CSR"){
        simDistsNN = lapply(integer(nSims), function(i){
            pSim = runifpoint(npTot, win = p$window)
            p1 = pSim[id <- sample(npTot, np1),];p2 = pSim[-id,]
            c(nncross(what = "dist", p1, p2), nncross(what = "dist", p2, p1))
        })
        mean(ecdf(unlist(simDistsNN))(obsDistNN))
    }
}
