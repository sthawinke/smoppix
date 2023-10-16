#' Estimate the probabilistic index for the collection of all distances
#'
#' @param pSub The subset point pattern wiht only one feature
#' @param ecdfAll The empirical cumulative distribution function of all distances under the null
#'
#' @return The probabilistic index for all distances
#' @importFrom spatstat.geom pairdist
calcAllDistPI = function(pSub, ecdfAll){
    obsDist = pairdist(pSub)
    mean(ecdfAll(obsDist))
}
