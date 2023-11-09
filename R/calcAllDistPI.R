#' Estimate the probabilistic index for the collection of all distances
#'
#' @param pSub The subset point pattern wiht only one feature
#' @param ecdfAll The empirical cumulative distribution function of all distances under the null
#'
#' @return The probabilistic index for all distances
#' @importFrom spatstat.geom pairdist
#' @inheritParams estPimsSingle
#' @importFrom stats dist
calcAllDistPI = function(pSub, p, ecdfAll, null, nSims){
    obsDist = dist(coords(pSub))
    if(null == "CSR"){
        mean(ecdfAll(obsDist))
    } else if(null == "background"){
        simDists = crossdist(pSub, subSampleP(p, nSims))
            # #Keep observed points fixed
        piEsts = vapply(seq_len(npoints(pSub)), FUN.VALUE = double(1), function(i){
            ecdf(simDists[i,])(obsDist[i])
        })
        mean(piEsts)
    }
}
calcAllDistPIpair = function(pSub1, pSub2, p, ecdfAll, null, nSims){
    obsDist = crossdist(pSub1, pSub2)
    if(null == "CSR"){
        mean(ecdfAll(obsDist))
    } else if(null == "background"){
        simDists1 = crossdist(pSub1, subSampleP(p, nSims))
        simDists2 = crossdist(pSub2, subSampleP(p, nSims))
        NP1 = npoints(pSub1);NP2 = npoints(pSub2)
        # #Keep observed points fixed
        piEsts1 = vapply(seq_len(NP1), FUN.VALUE = double(npoints(pSub2)), function(i){
            ecdf(simDists1[i,])(obsDist[i, ])
        })
        piEsts2 = vapply(seq_len(NP2), FUN.VALUE = double(npoints(pSub1)), function(i){
            ecdf(simDists2[i,])(obsDist[,i])
        })
        #Fix me
        mean(c(piEsts1, piEsts2))
    }
}
subSampleP = function(p, nSims){
    if(NP <- npoints(p) > nSims) p[sample(NP, nSims),] else p
}
