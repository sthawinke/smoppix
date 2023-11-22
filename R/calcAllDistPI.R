#' Estimate the probabilistic index for the collection of all distances
#'
#' @param pSub The subset point pattern with only one feature
#' @param ecdfAll The empirical cumulative distribution function of all distances under the null
#' @param rowSortMat A subset of distances, sorted along the rows
#'
#' @return The probabilistic index for all distances
#' @inheritParams estPimsSingle
#' @importFrom stats dist
calcAllDistPI = function(pSub, p, ecdfAll, null, rowSortMat){
    obsDist = as.matrix(dist(coords(pSub)))
    if(null == "CSR"){
        mean(ecdfAll(obsDist))
    } else if(null == "background"){
        piEsts = vapply(seq_len(nrow(rowSortMat)), FUN.VALUE = double(1), function(i){
            mean(ecdfPreSort(rowSortMat[i, ])(obsDist[i, -i]))
        })
        mean(piEsts)
    }
}
calcAllDistPIpair = function(NP1, NP2, ecdfAll, null, crossDist, rowSortMat){
    if(null == "CSR"){
        mean(ecdfAll(crossDist))
    } else if(null == "background"){
        # #Keep observed points fixed
        piEsts1 = sum(rowSortMat[seq_len(NP1), ] < crossDist)
        piEsts2 = sum(rowSortMat[-seq_len(NP1), ] < t(crossDist))
        (piEsts1+ piEsts2)/(NP1 +NP2)
    }
}

