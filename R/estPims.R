#' A wrapper function for the different probabilistic indices (PIs)
#'
#' @param p The point pattern
#' @param pis the Probabilitstic indices to be estimated
#' @param null A character vector, indicating how the null distribution is defined. See details.
#' @param nSims The number of Monte-Carlo simulations used to determine the null distribution if null = "CSR
#' @param nPointsAll How many points to subsample or simulate to calculate overall interpoint distance
#' and distance to point
#' @param nPointsAllWin How many points to subsample or simulate to calculate distance to cell edge or midpoint distribution
#' @param allowManyGenePairs A boolean, set to true to suppress warning messages for large numbers of gene pairs
#' @param manyPairs An integer, what are considered many gene pairs
#' @param verbose Should verbose output be printed?
#' @param owins The list of windows corresponding to cells
#' @param point The point to which the distances are to be calculated
#' @param features A character vector, for which features should the probablilistic indices be calculated?
#' @param tabObs A table of observed gene frequencies
#'
#' @return Data frames with estimated quantities per gene and/or gene pair
#' @importFrom stats ecdf dist
#' @importFrom spatstat.random runifpoint
#' @importFrom utils combn
#' @importFrom matrixStats colMins rowMins
#' @import spatstat.geom
#' @importFrom Rdpack reprompt
estPimsSingle = function(p, pis, null, tabObs, nSims = 5e1, nPointsAll = 2e3, nPointsAllWin = 2e2, features = NULL,
                   allowManyGenePairs = FALSE, manyPairs = 1e6, verbose = FALSE,
                   owins = NULL, point){
    if(is.null(features)){
        features = names(tabObs)
    } else {
        features = intersect(features, names(tabObs))
    }
    if(any(pis == "fixedpoint"))
        if(missing(point))
            stop("Provide fixed point to measure distances to when requesting pi = 'fixedPoint'")
        else
            pointPPP = ppp(point[1], point[2])
    if(any(pis %in% c("edge", "fixedpoint")) || (any(pis == "allDist") && null =="CSR")){
        if(any(pis %in% c("allDist", "fixedpoint"))){
           pSim = runifpoint(nPointsAll, win = p$window)
        }
        if(any(pis == "allDist"))
            ecdfAll = ecdf(dist(coords(pSim)))
        if(any(pis == "edge"))
            ecdfsEdge = lapply(owins, function(rr) ecdf(nncross(runifpoint(nPointsAllWin, win = rr), edges(rr), what = "dist")))
        if(any(pis == "midpoint"))
            ecdfMidPoint = lapply(owins, function(rr) ecdf(crossdist(runifpoint(nPointsAllWin, win = rr), centroid.owin(rr, as.ppp = TRUE))))
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
    uniPIs = lapply(nams <- names(tabObs[features])[!idOne], function(feat){
        pSub = p[id <- marks(p, drop = FALSE)$gene == feat, ]
        p = p[!id, ] #Avoid zero distances by removing observations of gene
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
    biPIs = if(any(piPair <- grepl(pis, pattern = "Pair"))){
        featuresBi = features[tabObs[features] > 0]
        if(!allowManyGenePairs && (numGenePairs <- choose(length(featuresBi), 2)) > manyPairs){
        warning(immediate. = TRUE, "Calculating probablistic indices for", numGenePairs, "gene pairs may take a long time!\n",
                "Set allowManyGenePairs to TRUE to suppress this message.")
        } else  if(verbose){
            message("Calculating bivariate probabilistic indices...")
        }
        genePairsMat = combn(featuresBi, 2)
        out = matrix(nrow = ncol(genePairsMat), byrow = TRUE, vapply(FUN.VALUE = double(sum(piPair)), seq_len(ncol(genePairsMat)), function(i){
            feat1 = genePairsMat[1, i];feat2 = genePairsMat[2, i]
            pSub1 = p[id1 <- which(marks(p, drop = FALSE)$gene == feat1), ]
            pSub2 = p[id2 <- which(marks(p, drop = FALSE)$gene == feat2), ]
            pJoin = p[c(id1, id2), ];p = p[-c(id1, id2), ]
            NNdistPI = if(any(pis == "nnPair") ){calcNNPIpair(pSub1, pSub2, p, pJoin = pJoin, null, nSims)}
            allDistPI = if(any(pis == "allDistPair")){
                calcAllDistPIpair(pSub1, pSub2, p, ecdfAll = ecdfAll, null = null, nSims = nPointsAll)}
            c("nnPair" = NNdistPI, "allDistPair" = allDistPI)
        }), dimnames = list(apply(genePairsMat, 2, paste, collapse = "--"), grep("Pair", pis, value = TRUE)))
        #Ensure it is a matrix even for one pi
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
#' @importFrom BiocParallel bplapply
#' @examples
#' data(Yang)
#' hypYang = suppressWarnings(buildHyperFrame(Yang, coordVars = c("x", "y"),
#' designVar = c("day", "root", "section")))
#' yangPims = estPims(hypYang, pis = c("nn", "nnPair"))
#' #Both univariate and bivariate tests, with nearest neighbour distances
#' @details
#' The null distribution used to calculate the PIs. Can be either "background",
#' in which case the observed distributions of all genes is used. Alternatively,
#' for null = "CSR", Monte-Carlo simulation under complete spatial randomness is performed within the given window.
#'
#' The 'nn' prefix indicates that nearest neighbour distances are being used, whereas 'all' indicates all distances are being used.
#' The suffix 'Pair' indicates that bivariate probabilistic indices, testing for co- and antilocalization are being used.
#' 'edge' and 'midpoint' calculate the distance to the edge respectively the midpoint of the windows added using the addCell() function.
#' 'fixedpoint' calculates the distances to a supplied list of points.
estPims = function(hypFrame, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint"),
                   null = c("background", "CSR"), features = attr(hypFrame, "features"),...){
    pis = match.arg(pis, several.ok = TRUE)
    null = match.arg(null)
    if(any(pis %in% c("edge", "midpoint")) && is.null(hypFrame$owins)){
        stop("No window provided for distance to edge or midpoint calculation. Add it using the addCell() function")
    }
    if(any(id <- !(features %in% attr(hypFrame, "features")))){
        stop("Features ", features[id], " not found in hyperframe")
    }
   out = bplapply(seq_len(nrow(hypFrame)), function(x){
       with(hypFrame[x,, drop = TRUE], estPimsSingle(ppp, owins = owins, pis = pis, null = null, tabObs = tabObs,...))
   })
   attr(out, "pis") = pis #Tag the pims calculated
   attr(out, "features") = features #Remember for which features the pims were calculated
   names(out) = rownames(hypFrame)
   out

}
