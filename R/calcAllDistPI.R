#' Estimate the probabilistic index for the collection of all distances
#'
#' @param pSub The subset point pattern wiht only one feature
#' @param ecdfAll The empirical cumulative distribution function of all distances under the null
#'
#' @return The probabilistic index for all distances
#' @importFrom spatstat.geom pairdist
#' @inheritParams estPims
calcAllDistPI = function(pSub, p, ecdfAll, null, nSims){
    obsDist = pairdist(pSub)
    if(null == "CSR"){
        mean(ecdfAll(obsDist))
    } else if(null == "background"){
        simDists = crossdist(pSub, if(npoints(p) > nSims) p[sample(npoints(p), nSims),] else p)
            # #Keep observed points fixed
        piEsts = vapply(seq_len(npoints(pSub)), FUN.VALUE = double(1), function(i){
            ecdf(simDists[i,])(obsDist[i])
        })
        mean(piEsts)
    }
}
