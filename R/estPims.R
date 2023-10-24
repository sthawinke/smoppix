#' A wrapper function for the different probabilistic indices (PIs)
#'
#' @param p The point pattern
#' @param pis the Probabilitstic indices to be estimated
#' @param ... Additional arguments passed on to the appropriate functions
#' @param null A character vector, indicating how the null distribution is defined. See details.
#' @param nSims The number of Monte-Carlo simulations used to determine the null distribution if null = "CSR
#' @param nPointsAll How many points to subsample to calculate overall distance distribution
#' @param allowManyGenePairs A boolean, set to true to suppress warning messages for large numbers of gene pairs
#' @param manyPairs An integer, what are considered many gene pairs
#' @param verbose Should verbose output be printed?
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
#' For the 'edge' and 'fixedpoint' probabilistic indices, always CSR is used for the null.
#'
#' @examples
estPimsSingle = function(p, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint"),
                   null = c("background", "CSR"), nSims = 1e2, nPointsAll = 1e4,
                   allowManyGenePairs = FALSE, manyPairs = 1e6, verbose = FALSE, rois = NULL, point,...){
    pis = match.arg(pis, several.ok = TRUE)
    if(any(pis == "edge") && is.null(rois)){
        stop("No region of interest provided for edge calculation. Add it using the addCell() function")
    }
    null = match.arg(null)
    tabObs = table(marks(p, drop = FALSE)$gene)
    unFeatures = names(tabObs); names(unFeatures) = unFeatures
    dfFeat = data.frame("feature" = unFeatures)
    if(any(pis %in% c("edge", "fixedpoint")) || (any(pis == "allDist") && null =="CSR")){
        if(any(pis == "allDist")){
            pSim = runifpoint(nPointsAll, win = p$window)
           ecdfAll = ecdf(dist(coords(pSim)))
        }
        if(any(pis == "edge"))
            ecdfsEdge = lapply(rois, function(rr) ecdf(nncross(runifpoint(nPointsAll, win = rr), rr, what = "dist")))
        if(any(pis == "fixedpoint"))
            ecdfPoint = ecdf(nncross(pSim, point, what = "dist"))
    }
    #Univariate patterns
    if(any(idZero <- (tabObs==1)) && verbose){
        message("Features\n", paste(sep = ", ", names(tabObs)[idZero]),
                "\nhave only one observations, so no intervent distances were calculated for them\n")
    }
    uniPIs = simplify2array(bplapply(unFeatures[!idZero], function(feat){
        pSub = subset(p, gene == feat)
        NNdistPI = if(any(pis == "nn") && (npoints(pSub) > 1)){calcNNPI(pSub, p, null, nSims)} else NULL
        allDistPI = if(any(pis == "allDist")&& (npoints(pSub) > 1)){
            calcAllDistPI(pSub, p, ecdfAll = ecdfAll, null = null, nSims = nPointsAll)} else NULL
        edgeDistPI = if(any(pis == "edge")){
            calcEdgeDistPI(pSub, p, ecdfAll = ecdfsEdge, ...)} else NULL
        pointDistPI = if(any(pis == "fixedpoint")){
            calcPointDistPI(pSub, p, ecdfAll = ecdfPoint, ...)} else NULL
        c("NNdistPI" = NNdistPI, "allDistPI" = allDistPI, "edgeDistPI" = edgeDistPI, "pointDistPI" = pointDistPI)
    }))
    #Bivariate patterns
    if(any(grepl(pis, pattern = "Pair"))){
        if(!allowManyGenePairs && (numGenePairs <- choose(length(unFeatures), 2)) > 1e6){
        warning(immediate. = TRUE, "Calculating probablistic indices for", numGenePairs, "gene pairs may take a long time!")
        }
    } else biPIs = NULL
    biPIs = NULL
    list("uniPIs" = uniPIs, "biPIs" = biPIs)
}
estPims = function(y, ...){
   with(y, estPimsSingle(ppp, rois = rois, ...))
}
