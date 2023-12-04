#' Calculate univariate PIs for all features
#'
#' @inheritParams estPimsSingle
#' @param ecdfsEdgeAndMidpoint ecdfs within the cells
#' @param ecdfAll Overall ecdf
#' @param ecdfs Ecdfs per point
#' @param idOne Index of features with only one observation, to be excluded
#' @param nSub The number of events in the subsampled pattern
#' @param verbose A boolean, should verbose output be printed
#' @param ecdfFixedPoint The ecdf of the fixed point
#' @param pointPPP the point pattern containing the fixed point
#'
#' @return PIs for every feature
#' @importFrom spatstat.geom marks
calcUniPIs = function(p, pis, verbose, ecdfsEdgeAndMidpoint, owins, tabObs, null, nSims,
                      ecdfs, nSub, ecdfAll, idOne, features, ecdfFixedPoint, pointPPP){
    if(verbose)
        message("Calculating univariate probabilistic indices...")
    uniPIs = lapply(nams <- names(tabObs[features])[!idOne], function(feat){
        pSub = p[id <-which(marks(p, drop = FALSE)$gene == feat), ]
        p = p[-id, ] #Avoid zero distances by removing observations of gene
        NNdistPI = if(any(pis == "nn")){
            calcNNPI(pSub, p, null, nSims, ecdfs = ecdfs[id], n = nSub,
                     ecdfAll = ecdfAll)}
        #Also here room for improvement
        allDistPI = if(any(pis == "allDist")){
            calcAllDistPI(pSub, p, ecdfAll = ecdfAll, null = null,
                          ecdfs = ecdfs[id])}
        edgeDistPI = if(any(pis == "edge")){
            calcWindowDistPI(pSub, owins, ecdfAll = ecdfsEdgeAndMidpoint,
                             towhat = "edge")}
        midPointDistPI = if(any(pis == "midpoint")){
            calcWindowDistPI(pSub, owins, ecdfAll = ecdfsEdgeAndMidpoint,
                             towhat = "midpoint")}
        nnCellPI = if(any(pis == "nnCell")){
            calcWindowDistPI(pSub, owins, ecdfAll = ecdfsEdgeAndMidpoint,
                             towhat = "nn")}
        allDistCellPI = if(any(pis == "allDistCell")){
            calcWindowDistPI(pSub, owins, ecdfAll = ecdfsEdgeAndMidpoint,
                             towhat = "allDist")}
        fixedPointDistPI = if(any(pis == "fixedpoint")){
            calcFixedPointDistPI(pSub, pointPPP, ecdfAll = ecdfFixedPoint)}
        list("pointDists" = c("nn" = NNdistPI, "allDist" = allDistPI),
             "windowDists" = list("edge" = edgeDistPI, "nnCell" = nnCellPI,
                                  "allDistCell" = allDistCellPI, "midpoint" = midPointDistPI,
                                  "fixedpoint" = fixedPointDistPI))
    }); names(uniPIs) = nams
    return(uniPIs)
}
