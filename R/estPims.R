#' A wrapper function for the different probabilistic indices (PIs)
#'
#' @param p The point pattern
#' @param pis The probabilistic indices to be estimated
#' @param null A character vector, indicating how the null distribution is
#'  defined. See details.
#' @param nPointsAll How many points to subsample or simulate to calculate
#' overall interpoint distance and distance to point.
#' This parameter is a strong driver of memory and cpu usage.
#' @param nPointsAllWin How many points to subsample or simulate
#' to calculate distance to cell edge or midpoint distribution
#' @param allowManyGenePairs A boolean, set to true to suppress
#' warning messages for large numbers of gene pairs
#' @param manyPairs An integer, what are considered many gene pairs
#' @param verbose Should verbose output be printed?
#' @param owins,centroids The list of windows corresponding to cells,
#' and their centroids
#' @param features A character vector, for which features should the
#' probabilistic indices be calculated?
#' @param tabObs A table of observed gene frequencies
#' @param verbose a boolean, should verbose output be printed?
#'
#' @return Data frames with estimated quantities per gene and/or gene pair
#' @importFrom stats ecdf dist
#' @importFrom spatstat.random runifpoint
#' @importFrom utils combn
#' @importFrom spatstat.geom nncross crossdist coords npoints
#' @importFrom Rdpack reprompt
#' @importFrom Rfast rowSort rowMins rowAny
#' @importFrom extraDistr pnhyper
estPimsSingle <- function(p, pis, null, tabObs, nPointsAll = 1e4,
    nPointsAllWin = 2e3, features = NULL, allowManyGenePairs = FALSE,
    manyPairs = 1e+06, verbose = FALSE, owins = NULL, centroids = NULL) {
    if (is.null(features)) {
        features <- names(tabObs)
    } else {
        features <- intersect(features, names(tabObs))
    }
    if (any(pis %in% c("nn", "nnPair"))) {
        if (null == "background") {
            # Subset some background points
            pSubLeft <- subSampleP(p, nPointsAll, returnId = TRUE)
        } else if (null == "CSR") {
            pSim <- runifpoint(nPointsAll, win = p$window)
            ecdfAll <- ecdf(stats::dist(coords(pSim)))
        }
    }
    if (any(pis %in% c("edge", "midpoint"))) {
        ecdfsCell <- lapply(names(owins), function(nam) {
            pSub <- switch(null,
                           "CSR" = runifpoint(nPointsAllWin, win = owins[[nam]]),
                           "background" = subSampleP(pCell <- p[marks(p, drop = FALSE)$cell == nam, ], nPointsAllWin)
            )
            edge <- if (any(pis == "edge")) {
                ecdf(nncross(pSub, edges(owins[[nam]]), what = "dist"))
            }
            midpoint <- if (any(pis == "midpoint")) {
                ecdf(crossdist(pSub, centroids[[nam]]))
            }
            list("edge" = edge, "midpoint" = midpoint)
        })
        names(ecdfsCell) <- names(owins)
    }
        NPall = npoints(p)
        pSplit = split.ppp(p, f = factor(marks(p, drop = FALSE)$gene))
        splitFac <- rep(seq_len(bpparam()$workers), length.out = length(nams <- names(tabObs[features])))
        #Divide the work over the available workers
        piList <- unsplit(f = splitFac, bplapply(split(nams, f = splitFac), function(ss){
            featPIs = lapply(ss, function(feat) {
                pSub <- pSplit[[feat]]
                NP <- npoints(pSub)
                distMat = cbind(
                    "self" = if(calcNNsingle <- (NP> 1 && ("nn" %in% pis))){
                        nndist(pSub)}, #Self distances
                        if("nnPair" %in% pis){
                            id = !(names(pSplit) %in% feat)
                    matrix(unlist(lapply(pSplit[id],
                        function(y){nncross(pSub, y, what = "dist")
                    })), nrow = NP, dimnames = list(NULL, names(pSplit)[id]))
                            #Cross-distances
                })
                if(is.matrix(distMat)){
                    approxRanksTmp = switch(null,
                        "background" = findRanksDist(getCoordsMat(pSub),
                            getCoordsMat(pSubLeft$Pout), distMat^2),
                        "CSR" = matrix(ecdfAll(distMat), nrow = nrow(distMat))
                    )
                    approxRanks = round(NPall*approxRanksTmp/(switch(null,
                        "background" = (npoints(pSubLeft$Pout)-
                             which(marks(p, drop = FALSE)$gene == feat) %in% pSubLeft$id),
                        #Correct for self distances by subtracting one.
                        #The C++ function only counts distances larger so self distances
                        #will be ignored in the numerator
                        "CSR" = 1)))
                    approxRanks[approxRanks==0] = 1
                    colnames(approxRanks) = colnames(distMat)
                }
                #Names may get lost in C++ function
                #And then rearrange to get to the PIs
                nnPI = if(calcNNsingle){
                     mean(calcNNPI(approxRanks[, "self"], NPall - (NP - 1),
                             m = NP - 1, r = 1))
                } else NA
                nnPIpair = if("nnPair" %in% pis){
                    apply(approxRanks[, colnames(approxRanks) != "self", drop = FALSE], 2,
                        calcNNPI, n = NPall - NP, m = NP, r = 1)
                }
                ## Window related distances
                edgeDistPI <- if (any(pis == "edge")) {
                    calcWindowDistPI(pSub, owins = owins, ecdfAll = ecdfsCell, pi = "edge")
                }
                midPointDistPI <- if (any(pis == "midpoint")) {
                    calcWindowDistPI(pSub, centroids = centroids,
                                     ecdfAll = ecdfsCell, pi = "midpoint")
                }
                list("windowDists" = list("edge" = edgeDistPI,
                                        "midpoint" = midPointDistPI),
                    "pointDists" = list("nn" = nnPI, "nnPair" = nnPIpair))
            })
            names(featPIs) = ss
            # TO DO: nnCell recycling code
            return(featPIs)
        }))
        names(piList) = nams
        nnPis = if("nn" %in% pis){
            vapply(piList[names(tabObs[tabObs>1])], FUN.VALUE = double(1), function(x){
                x$pointDists$nn
            })
        }
        genePairsMat <- combn(features, 2)
        if("nnPair" %in% pis){
            nnPairPis = vapply(seq_len(ncol(genePairsMat)), FUN.VALUE = double(1), function(i){
                feat1 <- genePairsMat[1, i]
                feat2 <- genePairsMat[2, i]
                mean(c(getElement(piList[[feat1]]$pointDists$nnPair, feat2),
                       getElement(piList[[feat2]]$pointDists$nnPair, feat1)))
            })
            names(nnPairPis) = apply(genePairsMat, 2, paste, collapse = "--")
        } else nnPairPis = NULL
        pointDists = list("nn" = nnPis, "nnPair" = nnPairPis)
        # Window pis
        windowDists = lapply(piList, function(x) x$windowDists)
        list("pointDists" = pointDists, "windowDists" = windowDists)
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
#' hypYang <- buildHyperFrame(Yang,
#'   coordVars = c('x', 'y'),
#'   imageVars = c('day', 'root', 'section')
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPims(hypYang[c(seq_len(5), seq(25, 29)),], pis = 'nn')
#' # Univariate nearest neighbour distances
#' @details
#' The null distribution used to calculate the PIs. Can be either 'background',
#' in which case the observed distributions of all genes is used. Alternatively,
#' for null = 'CSR', Monte-Carlo simulation under complete spatial randomness
#'  is performed within the given window.
#'
#' The 'nn' prefix indicates that nearest neighbour distances are being used,
#' whereas 'all' indicates all distances are being used.
#' The suffix 'Pair' indicates that bivariate probabilistic indices,
#' testing for co- and antilocalization are being used.
#' 'edge' and 'midpoint' calculate the distance to the edge respectively
#'  the midpoint of the windows added using the addCell() function.
#' The suffix 'Cell' indicates distances are being calculated within cells only.
estPims <- function(hypFrame, pis = c("nn", "nnPair",
    "edge", "midpoint", "nnCell", "nnPairCell"), verbose = TRUE,
    null = c("background", "CSR"),
    features = getFeatures(hypFrame),...) {
    pis <- match.arg(pis, several.ok = TRUE)
    null <- match.arg(null)
    if (any(pis %in% c("edge", "midpoint", "nnCell")) &&
        is.null(hypFrame$owins)) {
        stop("No window provided for distance to edge or midpoint calculation. "
             , "Add it using the addCell() function")
    }
    if (any(id <- !(features %in% getFeatures(hypFrame)))) {
        stop("Features ", features[id], " not found in hyperframe")
    }
    if (verbose)
        message("Calculating PIs for point pattern ")
    hypFrame$pimRes <- lapply(seq_len(nrow(hypFrame)), function(x) {
        if (verbose) {
            message(x, " of ", nrow(hypFrame), "  ")
        }
       out <- estPimsSingle(hypFrame[[x, "ppp"]],
                      owins = hypFrame[x, "owins", drop = TRUE], pis = pis,
                      null = null, tabObs = hypFrame[[x,"tabObs"]],
                      centroids = hypFrame[x, "centroids", drop = TRUE], ...)
       return(out)
    })
    list(hypFrame = hypFrame, null = null, pis = pis, features = features)
}
