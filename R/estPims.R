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
#' @import spatstat.geom
#' @importFrom Rdpack reprompt
#' @importFrom Rfast rowSort rowMins rowAny
#' @importFrom extraDistr pnhyper
estPimsSingle = function(p, pis, null, tabObs, nSims = 5e1, nPointsAll = 2e3,
                         nPointsAllWin = 5e2, features = NULL,
                   allowManyGenePairs = FALSE, manyPairs = 1e6, verbose = FALSE,
                   owins = NULL, point){
    if(is.null(features)){
        features = names(tabObs)
    } else {
        features = intersect(features, names(tabObs))
    }
    if(any(pis == "fixedpoint"))
        if(missing(point))
            stop("Provide fixed point to measure distances to when requesting",
                 "pi = 'fixedPoint'")
        else
            pointPPP = ppp(point[1], point[2])
    if(any(pis %in% c("edge", "fixedpoint", "allDist")) && null =="CSR"){
        if(any(pis %in% c("allDist", "fixedpoint"))){
           pSim = runifpoint(nPointsAll, win = p$window)
        }
        if(any(pis == "allDist"))
            ecdfAll = ecdf(dist(coords(pSim)))
        if(any(pis %in% c("edge", "midpoint", "nnCell", "allDistCell"))){
            ecdfsEdgeAndMidpoint = lapply(owins, function(rr) {
                pSub = runifpoint(nPointsAllWin, win = rr)
                edge = if(any(pis == "edge"))
                    ecdf(nncross(pSub, edges(rr), what = "dist"))
                midpoint = if(any(pis == "midpoint"))
                    ecdf(crossdist(pSub, centroid.owin(rr, as.ppp = TRUE)))
                nnCell = if(any(pis == "nnCell"))
                    ecdf(nndist(pSub))
                allDistCell = if(any(pis == "allDistCell"))
                    ecdf(dist(coords(pSub)))
                list("edge" = edge, "midpoint" = midpoint, "nnCell" = nnCell,
                     "allDistCell" = allDistCell)
            });names(ecdfsEdgeAndMidpoint) = names(owins)
        }
        if(any(pis == "fixedpoint"))
            ecdfFixedPoint = ecdf(nncross(pSim, pointPPP, what = "dist"))
    } else if(any(pis %in% c("edge", "fixedpoint", "midpoint")) &&
              null =="background"){
        if(any(pis %in% c("edge", "midpoint"))){
            ecdfsEdgeAndMidpoint = lapply(names(owins), function(nam) {
                pSub = subSampleP(p[marks(p, drop = FALSE)$cell == nam,],
                                  nPointsAllWin)
                edge = if(any(pis == "edge"))
                    ecdf(nncross(pSub, edges(owins[[nam]]), what = "dist"))
                midpoint = if(any(pis == "midpoint"))
                    ecdf(crossdist(pSub,
                                   centroid.owin(owins[[nam]], as.ppp = TRUE)))
                list("edge" = edge, "midpoint" = midpoint)
                });names(ecdfsEdgeAndMidpoint) = names(owins)
        }
        if(any(pis == "fixedpoint"))
            ecdfFixedPoint = ecdf(nncross(subSampleP(p, nPointsAll),
                                          pointPPP, what = "dist"))
    }
    if(any(pis %in% c("nn", "nnPair", "allDist", "allDistPair"))){
            #Prepare some null distances, either with CSR or background
            subSamSize = max(c(nPointsAll, tabObs+1e2))
            pSub = switch(null,
                          "background" = subSampleP(p, subSamSize),
                          "CSR" = runifpoint(subSamSize, win = p$window))
            #Subsample for memory reasons
            cd = getCrossDist(p, pSub, null)
            ecdfs = apply(rowSort(cd), 1, ecdfPreSort)
            nSub = ncol(cd)
    }
    #Univariate patterns
    if(any(idOne <- (tabObs[features]==1)) && verbose){
        message("Features\n", paste(collapse = ", ", names(tabObs[features])[idOne]),
                "\nhave only one observation, so no intervent distances were calculated for them\n")
    }
    piPair <- grepl(pis, pattern = "Pair")
    uniPIs = if(any(!piPair)){
        calcUniPIs(p, pis, verbose, ecdfsEdgeAndMidpoint, owins, tabObs, null,
                   nSims, ecdfs, nSub, ecdfAll, idOne, features, ecdfFixedPoint,
                   pointPPP)
    }
    #Bivariate patterns
    biPIs = if(any(piPair)){
        calcBiPIs(features = features[tabObs[features] > 0],p, pis, null, nSims,
                  ecdfs, nSub, ecdfAll, manyPairs, verbose, allowManyGenePairs)
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
#' hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
#' imageVars = c("day", "root", "section"))
#' #Fit a subset of features to limit computation time
#' yangPims = estPims(hypYang, pis = c("nn", "nnPair"))
#' #Both univariate and bivariate tests, with nearest neighbour distances
#' @details
#' The null distribution used to calculate the PIs. Can be either "background",
#' in which case the observed distributions of all genes is used. Alternatively,
#' for null = "CSR", Monte-Carlo simulation under complete spatial randomness
#'  is performed within the given window.
#'
#' The 'nn' prefix indicates that nearest neighbour distances are being used,
#' whereas 'all' indicates all distances are being used.
#' The suffix 'Pair' indicates that bivariate probabilistic indices,
#' testing for co- and antilocalization are being used.
#' 'edge' and 'midpoint' calculate the distance to the edge respectively
#'  the midpoint of the windows added using the addCell() function.
#' 'fixedpoint' calculates the distances to a supplied list of points.
#' The suffix "Cell" indicates distances are being calculated within cells only.
estPims = function(hypFrame, pis = c("nn", "allDist", "nnPair", "allDistPair",
                    "edge", "midpoint", "fixedpoint", "nnCell", "allDistCell"),
                   null = c("background", "CSR"),
                   features = attr(hypFrame, "features"),...){
    pis = match.arg(pis, several.ok = TRUE)
    null = match.arg(null)
    if(any(pis %in% c("edge", "midpoint", "nnCell", "allDistCell")) && is.null(hypFrame$owins)){
        stop("No window provided for distance to edge or midpoint calculation.",
             "Add it using the addCell() function")
    }
    if(any(id <- !(features %in% attr(hypFrame, "features")))){
        stop("Features ", features[id], " not found in hyperframe")
    }
    hypFrame$pimRes = bplapply(seq_len(nrow(hypFrame)), function(x){
       with(hypFrame[x,, drop = TRUE], estPimsSingle(ppp, owins = owins,
            pis = pis, null = null, tabObs = tabObs,...))
   })
   attr(hypFrame, "pis") = pis #Tag the pims calculated
   attr(hypFrame, "featuresEst") = features
   #Remember for which features the pims were calculated
   hypFrame
}
