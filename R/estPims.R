#' A wrapper function for the different probabilistic indices (PIs)
#'
#' @param p The point pattern
#' @param pis the Probabilitstic indices to be estimated
#' @param null A character vector, indicating how the null distribution is defined. See details.
#' @param nSims The number of Monte-Carlo simulations used to determine the null distribution if null = "CSR
#' @param nPointsAll How many points to subsample or simulate to calculate overall interpoint distance
#' and distance to point or edge distribution
#' @param allowManyGenePairs A boolean, set to true to suppress warning messages for large numbers of gene pairs
#' @param manyPairs An integer, what are considered many gene pairs
#' @param verbose Should verbose output be printed?
#' @param owins The list of windows corresponding to cells
#' @param point The point to which the distances are to be calculated
#' @param features A character vector, for which features should the probablilistic indices be calculated?
#' @param tabObs A table of observed gene frequencies
#'
#' @return Data frames with estimated quantities per gene and/or gene pair
#' @export
#' @importFrom BiocParallel bplapply
#' @importFrom stats ecdf dist
#' @importFrom spatstat.random runifpoint
#' @importFrom utils combn
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
estPimsSingle = function(p, pis, null, tabObs, nSims = 5e1, nPointsAll = 2e3, features = NULL,
                   allowManyGenePairs = FALSE, manyPairs = 1e6, verbose = FALSE,
                   owins = NULL, point){
    if(is.null(features))
        features = names(tabObs)
    else
        features = intersect(features, names(tabObs))
    if(any(pis == "fixedpoint"))
        pointPPP = ppp(point[1], point[2])
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
            ecdfFixedPoint = ecdf(nncross(pSim, pointPPP, what = "dist"))
    } else if(any(pis %in% c("edge", "fixedpoint", "midpoint")) && null =="background"){
        pSubAll = subSampleP(p, nPointsAll)
        if(any(pis == "edge"))
            ecdfsEdge = lapply(owins, function(rr) ecdf(nncross(pSubAll, edges(rr), what = "dist")))
        if(any(pis == "midpoint"))
            ecdfMidPoint = lapply(owins, function(rr) ecdf(crossdist(pSubAll, centroid.owin(rr, as.ppp = TRUE))))
        if(any(pis == "fixedpoint"))
            ecdfFixedPoint = ecdf(nncross(pSubAll, pointPPP, what = "dist"))
    }
    #Univariate patterns
    if(any(idOne <- (tabObs[features]==1)) && verbose){
        message("Features\n", paste(sep = ", ", names(tabObs[features])[idOne]),
                "\nhave only one observation, so no intervent distances were calculated for them\n")
    }
    if(verbose)
        message("Calculating univariate probabilistic indices...")
    uniPIs = bplapply(nams <- names(tabObs[features])[!idOne], function(feat){
        pSub = subset(p, gene == feat)
        NNdistPI = if(any(pis == "nn") ){calcNNPI(pSub, p, null, nSims)}
        allDistPI = if(any(pis == "allDist")){
            calcAllDistPI(pSub, p, ecdfAll = ecdfAll, null = null, nSims = nPointsAll)}
        edgeDistPI = if(any(pis == "edge")){
            calcWindowDistPI(pSub, owins, ecdfAll = ecdfsEdge, towhat = "edge")}
        midPointDistPI = if(any(pis == "midpoint")){
            calcWindowDistPI(pSub, owins, ecdfAll = ecdfMidPoint, towhat = "midpoint")}
        fixedPointDistPI = if(any(pis == "fixedpoint")){
            calcFixedPointDistPI(pSub, pointPPP, ecdfAll = ecdfFixedPoint)}
        list("pointDists" = c("nn" = NNdistPI, "allDist" = allDistPI),
             "windowDists" = list("edge" = edgeDistPI, "midpoint" = midPointDistPI, "fixedpoint" = fixedPointDistPI))
    }); names(uniPIs) = nams
    #Bivariate patterns
    biPIs = if(any(grepl(pis, pattern = "Pair"))){
        featuresBi = features[tabObs[features] > 0]
        if(!allowManyGenePairs && (numGenePairs <- choose(length(featuresBi), 2)) > 1e6){
        warning(immediate. = TRUE, "Calculating probablistic indices for", numGenePairs, "gene pairs may take a long time!\n",
                "Set allowManyGenePairs to TRUE to suppress this message.")
        } else  if(verbose){
            message("Calculating bivariate probabilistic indices...")
        }
        genePairsMat = combn(featuresBi, 2)
        out = simplify2array(bplapply(seq_len(ncol(genePairsMat)), function(i){
            feat1 = genePairsMat[1, i];feat2 = genePairsMat[2, i]
            pSub1 = subset(p, gene == feat1)
            pSub2 = subset(p, gene == feat2)
            p = subset(p, !(gene %in% genePairsMat[, i]))
            NNdistPI = if(any(pis == "nn") ){calcNNPIpair(pSub1, pSub2, p, null, nSims)}
            allDistPI = if(any(pis == "allDist")){
                calcAllDistPIpair(pSub1, pSub2, p, ecdfAll = ecdfAll, null = null, nSims = nPointsAll)}
           # NPsort = sort(tabObs[genePairsMat[, i]]); names(NPsort) = c("minNP", "maxNP")
            c("nnPair" = NNdistPI, "allDistPair" = allDistPI)
        }))
        colnames(out) = apply(genePairsMat, 2, paste, collapse = "_")
        out
    }
    list("uniPIs" = uniPIs, "biPIs" = biPIs)
}
#' Estimate pims on a hyperframe
#' @export
#' @param hypFrame the hyperframe
#' @param ... additional arguments, passed on to estPimsSingle
#' @inheritParams estPimsSingle
#' @return A list of estimated pims
#' @examples
estPims = function(hypFrame, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint"),
                   null = c("background", "CSR"), ...){
    pis = match.arg(pis, several.ok = TRUE)
    null = match.arg(null)
    if(any(pis %in% c("edge", "midpoint")) && is.null(hypFrame$owins)){
        stop("No window provided for distance to edge or midpoint calculation. Add it using the addCell() function")
    }
   out = with(hypFrame, estPimsSingle(ppp, owins = owins, pis = pis, null = null, tabObs = table,...))
   attr(out, "pis") = pis
   #attr(out, "features") = features
   out

}
