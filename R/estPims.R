#' A wrapper function for the different probabilistic indices (PIs)
#'
#' @param p The point patterns
#' @param pis the Probabilitstic indices to be estimated
#' @param ... Additional arguments passed on to the appropriate functions
#' @param null A character vector, indicating how the null distribution is defined. See details.
#' @param nSims The number of Monte-Carlo simulations used to determine the null distribution if null = "CSR
#' @param nPointsAll How many points to subsample to calculate overall distance distribution
#' @param allowManyGenePairs A boolean, set to true to suppress warning messages for large numbers of gene pairs
#' @param manyPairs An integer, what are considered many gene pairs
#'
#' @return Data frames with estimated quantities per gene and/or gene pair
#' @export
#' @importFrom BiocParallel bplapply
#' @importFrom stats ecdf dist
#' @import spatstat.geom
#' @importFrom spatstat.random runifpoint
#' @details
#' The null distribution used to calculate the PIs. Can be either "background",
#' in which case the observed distributions of all genes is used. Alternatively,
#' for null = "CSR", Monte-Carlo simulation under complete spatial randomness is performed within the given window.
#'
#' @examples
estPims = function(p, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint"),
                   null = c("background", "CSR"), nSims = 1e2, nPointsAll = 1e4,
                   allowManyGenePairs = FALSE, manyPairs = 1e6, ...){
    pis = match.arg(pis, several.ok = TRUE)
    null = match.arg(null)
    tabObs = table(marks(p)$gene)
    unFeatures = names(tabObs); names(unFeatures) = unFeatures
    dfFeat = data.frame("feature" = unFeatures)
    if(any(pis == "allDist") && null =="CSR"){
        ecdfAll = ecdf(dist(coords(runifpoint(nPointsAll, win = p$window))))
        #"background" = if(np <- npoints(p) > nPointsAll) p[sample(np, nPointsAll)] else p)
    }
    #Univariate patterns
    uniPIs = bplapply(unFeatures, function(feat){
        pSub = subset(p, gene == feat)
        NNdistPI = if(any(pis == "nn")){calcNNPI(pSub, p, null, nSims)} else NULL
        allDistPI = if(any(pis == "allDist")){calcAllDistPI(pSub, p, ecdfAll = ecdfAll, null = null, nSims = nPointsAll)} else NULL
        edgeDistPI = if(any(pis == "edge")){calcEdgeDistPI(pSub, null, nSims, ...)} else NULL
        pointDistPI = if(any(pis == "fixedpoint")){calcPointDistPI(pSub, null, nSims, ...)} else NULL
        c("NNdistPI" = NNdistPI, "allDistPI" = allDistPI, "edgeDistPI" = edgeDistPI, "pointDistPI" = pointDistPI)
    })
    #Bivariate patterns
    if(any(grepl(pis, "genePair")) && !allowManyGenePairs && (numGenePairs <- choose(length(unFeatures), 2)) > 1e6){
        warning(immediate. = TRUE, "Calculating probablistic indices for", numGenePairs, "gene pairs may take a long time!")
    }


}
