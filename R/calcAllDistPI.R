#' Estimate the probabilistic index for the collection of all distances
#'
#' @param pSub The subset point pattern with only one feature
#' @param ecdfAll The empirical cumulative distribution function of all distances under the null
#'
#' @return The probabilistic index for all distances
#' @inheritParams estPimsSingle
#' @importFrom stats dist
calcAllDistPI = function(pSub, p, ecdfAll, null, nSims, id, rowSortMat){
    obsDist = as.matrix(dist(coords(pSub)))
    if(null == "CSR"){
        mean(ecdfAll(obsDist))
    } else if(null == "background"){
        # simDists = crossdist(pSub, subSampleP(p, nSims))
        #     # #Keep observed points fixed
        # piEsts = vapply(seq_len(npoints(pSub)), FUN.VALUE = double(1), function(i){
        #     mean(simDists[i,] < obsDist[i,-i])
        # })
        # mean(piEsts)
        piEsts = vapply(seq_along(id), FUN.VALUE = double(1), function(i){
            mean(ecdfPreSort(rowSortMat[id[i], ])(obsDist[i, -i]))
        })
        mean(piEsts)
    }
}
calcAllDistPIpair = function(NP1, NP2, ecdfAll, null, crossDist, backDist){
    if(null == "CSR"){
        mean(ecdfAll(crossDist))
    } else if(null == "background"){
        # #Keep observed points fixed
        piEsts1 = sum(backDist[seq_len(NP1), ] < crossDist)
        piEsts2 = sum(backDist[-seq_len(NP1), ] < t(crossDist))
        (piEsts1+ piEsts2)/(NP1 +NP2)
    }
}

