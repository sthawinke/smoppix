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
    if(verbose)
        message("Calculating univariate probabilistic indices...")
    uniPIs = lapply(nams <- names(tabObs[features])[!idOne], function(feat){
        pSub = p[id <-which(marks(p, drop = FALSE)$gene == feat), ]
        p = p[-id, ] #Avoid zero distances by removing observations of gene
        NNdistPI = if(any(pis == "nn")){
            calcNNPI(pSub, p, null, nSims, ecdfs = ecdfs[id], n = nSub)}
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
    #Bivariate patterns
    biPIs = if(any(piPair <- grepl(pis, pattern = "Pair"))){
        featuresBi = features[tabObs[features] > 0]
        if(!allowManyGenePairs &&
           (numGenePairs <- choose(length(featuresBi), 2)) > manyPairs){
        warning(immediate. = TRUE, "Calculating probablistic indices for",
                numGenePairs, "gene pairs may take a long time!\n",
                "Set allowManyGenePairs to TRUE to suppress this announcement")
        } else  if(verbose){
            message("Calculating bivariate probabilistic indices...")
        }
        genePairsMat = combn(featuresBi, 2)
        out = matrix(nrow = ncol(genePairsMat), byrow = TRUE,
            vapply(seq_len(ncol(genePairsMat)), FUN.VALUE = double(sum(piPair)),
                   function(i){
            feat1 = genePairsMat[1, i];feat2 = genePairsMat[2, i]
            pSub1 = p[id1 <- which(marks(p, drop = FALSE)$gene == feat1), ]
            pSub2 = p[id2 <- which(marks(p, drop = FALSE)$gene == feat2), ]
            cd = crossdist(pSub1, pSub2)
            #Reorder and subset if needed
            NNdistPI = if(any(pis == "nnPair")){
                calcNNPIpair(cd = cd, id1 = id1, id2 = id2, null = null,
                             ecdfs = ecdfs[c(id1, id2)], p = p,
                             nSims = nSims, n = nSub)
                }
            allDistPI = if(any(pis == "allDistPair")){
                calcAllDistPIpair(id1 = id1, id2 = id2, ecdfAll = ecdfAll,
                                  null = null, ecdfs = ecdfs[c(id1, id2)],
                                  crossDist = cd)
                }
            c("nnPair" = NNdistPI, "allDistPair" = allDistPI)
        }), dimnames = list(apply(genePairsMat, 2, paste, collapse = "--"),
                            grep("Pair", pis, value = TRUE)))
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
