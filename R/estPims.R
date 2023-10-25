#' A wrapper function for the different probabilistic indices (PIs)
#'
#' @param p The point pattern
#' @param pis the Probabilitstic indices to be estimated
#' @param ... Additional arguments passed on to the appropriate functions
#' @param null A character vector, indicating how the null distribution is defined. See details.
#' @param nSims The number of Monte-Carlo simulations used to determine the null distribution if null = "CSR
#' @param nPointsAll How many points to subsample or simulate to calculate overall interpoint distance
#' and distance to point or edge distribution
#' @param allowManyGenePairs A boolean, set to true to suppress warning messages for large numbers of gene pairs
#' @param manyPairs An integer, what are considered many gene pairs
#' @param verbose Should verbose output be printed?
#' @param points A list of points to which the distances are to be calculated
#'
#' @return Data frames with estimated quantities per gene and/or gene pair
#' @export
#' @importFrom BiocParallel bplapply
#' @importFrom stats ecdf dist
#' @importFrom spatstat.random runifpoint
#' @import spatstat.geom
#' @details
#' The null distribution used to calculate the PIs. Can be either "background",
#' in which case the observed distributions of all genes is used. Alternatively,
#' for null = "CSR", Monte-Carlo simulation under complete spatial randomness is performed within the given window.
#' For the 'edge', 'midpoint' and 'fixedpoint' probabilistic indices, always CSR is used for the null.
#'
#' The 'nn' prefix indicates that nearest neighbour distances are being used, whereas 'all' indicates all distances are being used.
#' The suffix 'Pair' indicates that bivariate probabilistic indices, testing for co- and antilocalization are being used.
#' 'edge' and 'midpoint' calculate the distance to the edge respectively the midpoint of the windows added using the addCell() function.
#' 'fixedpoint' calculates the distances to a supplied list of points.
#' @examples
estPimsSingle = function(p, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint"),
                   null = c("background", "CSR"), nSims = 1e2, nPointsAll = 5e3,
                   allowManyGenePairs = FALSE, manyPairs = 1e6, verbose = FALSE, owins = NULL, point,...){
    pis = match.arg(pis, several.ok = TRUE)
    if(any(pis == "edge") && is.null(owins)){
        stop("No region of interest provided for edge calculation. Add it using the addCell() function")
    }
    if(any(pis %in% c("edge", "fixedpoint", "midpoint")) && null == "background"){
        message("For distance to edge, cell centroid and fixed point, only CSR is implemented as null and will be used")
    }
    null = match.arg(null)
    tabObs = table(marks(p, drop = FALSE)$gene)
    unFeatures = names(tabObs); names(unFeatures) = unFeatures
    dfFeat = data.frame("feature" = unFeatures)
    if(any(pis %in% c("edge", "fixedpoint")) || (any(pis == "allDist") && null =="CSR")){
        if(any(pis %in% c("allDist", "fixedpoint"))){
           pSim = runifpoint(nPointsAll, win = p$window)
        }
        if(any(pis == "allDist"))
            ecdfAll = ecdf(dist(coords(pSim)))
        if(any(pis == "edge"))
            ecdfsEdge = lapply(owins, function(rr) ecdf(nncross(runifpoint(nPointsAll, win = rr), edges(rr), what = "dist")))
        if(any(pis == "midpoint"))
            ecdfMidPoint = lapply(owins, function(rr) ecdf(crossdist(runifpoint(nPointsAll, win = rr), centroid.owin(rr, as.ppp = TRUE))))
        if(any(pis == "fixedpoint"))
            pointPPP = ppp(point[1], point[2])
            ecdfFixedPoint = ecdf(nncross(pSim, pointPPP, what = "dist"))
    }
    #Univariate patterns
    if(any(idZero <- (tabObs==1)) && verbose){
        message("Features\n", paste(sep = ", ", names(tabObs)[idZero]),
                "\nhave only one observations, so no intervent distances were calculated for them\n")
    }
    uniPIs = bplapply(unFeatures[!idZero], function(feat){
        pSub = subset(p, gene == feat)
        NNdistPI = if(any(pis == "nn") && (npoints(pSub) > 1)){calcNNPI(pSub, p, null, nSims)}
        allDistPI = if(any(pis == "allDist")&& (npoints(pSub) > 1)){
            calcAllDistPI(pSub, p, ecdfAll = ecdfAll, null = null, nSims = nPointsAll)}
        edgeDistPI = if(any(pis == "edge")){
            calcWindowDistPI(pSub, owins, ecdfAll = ecdfsEdge)}
        midPointDistPI = if(any(pis == "midpoint")){
            calcWindowDistPI(pSub, owins, ecdfAll = ecdfMidPoint, midPoint = TRUE)}
        fixedPointDistPI = if(any(pis == "fixedpoint")){
            calcFixedPointDistPI(pSub, pointPPP, ecdfAll = ecdfFixedPoint)}
        list("pointDists" = c("NNdistPI" = NNdistPI, "allDistPI" = allDistPI, "fixedPointDistPI" = fixedPointDistPI),
             "windowDists" = list("edgeDistPI" = edgeDistPI, "midPointDistPI" = midPointDistPI))
    })
    #Bivariate patterns
    if(any(grepl(pis, pattern = "Pair"))){
        if(!allowManyGenePairs && (numGenePairs <- choose(length(unFeatures), 2)) > 1e6){
        warning(immediate. = TRUE, "Calculating probablistic indices for", numGenePairs, "gene pairs may take a long time!\n",
                "Set allowManyGenePairs to TRUE to suppress this message.")
        }
    } else biPIs = NULL
    biPIs = NULL
    list("uniPIs" = uniPIs, "biPIs" = biPIs)
}
#' Estimate pims on a hyperframe
#' @export
#' @param hypFrame the hyperframe
#' @return A list of estimated pims
#' @examples
estPims = function(hypFrame, ...){
   with(hypFrame, estPimsSingle(ppp, owins = owins, ...))
}
