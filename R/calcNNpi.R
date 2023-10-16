#' Estimate the PI for the nearest neighbour distances
#'
#' @param pSub The subset point pattern containing only a single gene
#' @inheritParams estPims
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
            #nndist(p[sample(npoints(p), npoints(p)),])
            distMat = crossdist(pSub, p[sample(npoints(p), npoints(p)-1),])
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
        simDistsNN = vapply(integer(nSims), FUN.VALUE = double(npoints(pSub)), function(i){
            nndist(rpoispp(lambda, win = pSub$window))
        })
        mean(ecdf(simDistsNN)(obsDistNN))
    }
}
