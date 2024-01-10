#' A wrapper function for the different probabilistic indices (PIs)
#'
#' @param p The point pattern
#' @param pis The probabilistic indices to be estimated
#' @param null A character vector, indicating how the null distribution is
#'  defined. See details.
#' @param nPointsAll,nPointsAllWithinCell How many points to subsample or simulate to calculate
#' overall interpoint distance and distance to point.
#' This parameter is a strong driver of computation time.
#'The second argument applies to within cell calculations.
#' @param nPointsAllWin How many points to subsample or simulate
#' to calculate distance to cell edge or midpoint distribution
#' @param owins,centroids The list of windows corresponding to cells,
#' and their centroids
#' @param features A character vector, for which features should the
#' probabilistic indices be calculated?
#' @param tabObs A table of observed gene frequencies
#' @param window An window of class owin, in which events can occur
#' @param loopFun The function to use to loop over the features.
#' Defaults to bplapply except when looping over features within cells
#' @param minDiff An integer, the minimum number of events from other genes
#'  needed for calculation of background PIs
#'
#' @return Data frames with estimated quantities per gene and/or gene pair
#' @importFrom stats ecdf dist
#' @importFrom spatstat.random runifpoint
#' @importFrom Rdpack reprompt
#' @importFrom BiocParallel bplapply
estPimsSingle <- function(p, pis, null, tabObs, nPointsAll = switch(null, background = 1e5, CSR = 1e3),
                          nPointsAllWithinCell = switch(null, background = 1e4, CSR = 5e2), nPointsAllWin = 1e3,
    features = NULL, owins = NULL, centroids = NULL, window = p$window, loopFun = "bplapply", minDiff = 2e1) {
    if(!length(features)){
        return(list(pointDists = NULL, windowDists = NULL, withinCellDists = NULL))
    }
    features <- sample(intersect(features, names(tabObs)))
    #Scramble to ensure equal calculation times in multithreading
    if (any(pis %in% c("nn", "nnPair"))) {
        if (null == "background") {
            # Subset some background points
            pSubLeft <- subSampleP(p, nPointsAll, returnId = TRUE)
        } else if (null == "CSR") {
            pSim <- runifpoint(nPointsAll, win = window)
            ecdfAll <- ecdf(stats::dist(coords(pSim)))
        }
    }
    if (any(pis %in% c("edge", "midpoint"))) {
        ecdfsCell <- findEcdfsCell(p = p, owins = owins, centroids = centroids,
                nPointsAllWin = nPointsAllWin, null = null, pis = pis)
    }
    piList <- calcIndividualPIs(p = p, pSubLeft = pSubLeft, pis = pis,
                                null = null, tabObs = tabObs, owins = owins,
        centroids = centroids, features = features, ecdfAll = ecdfAll,
        ecdfsCell = ecdfsCell, loopFun = loopFun, minDiff = minDiff)
    nnPis <- if ("nn" %in% pis) {
        vapply(piList[intersect(features, names(tabObs[tabObs > 1]))],
               FUN.VALUE = double(1), function(x) {
            x$pointDists$nn
        })
    }
    if ("nnPair" %in% pis && length(features) > 1 ) {
        genePairsMat <- combn(features, 2)
        nnPairPis <- vapply(seq_len(ncol(genePairsMat)), FUN.VALUE = double(1),
                            function(i) {
            feat1 <- genePairsMat[1, i]
            feat2 <- genePairsMat[2, i]
            mean(c(getElement(piList[[feat1]]$pointDists$nnPair, feat2),
                   getElement(piList[[feat2]]$pointDists$nnPair,
                feat1)))
        })
        names(nnPairPis) <- apply(genePairsMat, 2, paste, collapse = "--")
    } else nnPairPis <- NULL
    pointDists <- list(nn = nnPis, nnPair = nnPairPis)
    # Window pis
    windowDists <- if (any(pis %in% c("edge", "midpoint"))) {
        lapply(piList, function(x) x$windowDists)
    }
    # Within cell: recurse into estPimsSingle but now per cell
    withinCellDists <- if (any(grepl("Cell", pis))) {
        unCells <- setdiff(unique(marks(p, drop = FALSE)$cell), "NA")
        cellDists <- bplapply(unCells, function(nam) {
            pSub = p[marks(p, drop = FALSE)$cell == nam,]
            estPimsSingle(pis = gsub("Cell", "", grep(value = TRUE, "Cell", pis)),
                          p = pSub, null = null,
                nPointsAll = nPointsAllWithinCell, window = owins[[nam]], features = features,
                tabObs = table(marks(pSub, drop = FALSE)$gene),
                loopFun = "lapply")$pointDists
        })
        names(cellDists) <- unCells
        cellDists
    }
    list(pointDists = pointDists, windowDists = windowDists,
         withinCellDists = withinCellDists)
}
#' Estimate probabilistic indices for (co)localization on a hyperframe
#' @description Estimate different probabilistic indices for (co)localization
#' on  a hyperframe, and integrate the results in the same hyperframe
#' @export
#' @param hypFrame the hyperframe
#' @param ... additional arguments, passed on to estPimsSingle
#' @param verbose a boolean, whether to report on progress of fitting
#' @inheritParams estPimsSingle
#' @return The hyperframe with the estimated PIMs present in it
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
estPims <- function(hypFrame, pis = c("nn", "nnPair", "edge", "midpoint", "nnCell", "nnPairCell"), verbose = TRUE, null = c("background",
    "CSR"), features = getFeatures(hypFrame), ...) {
    pis <- match.arg(pis, several.ok = TRUE)
    null <- match.arg(null)
    if (any(pis %in% c("edge", "midpoint", "nnCell")) && is.null(hypFrame$owins)) {
        stop("No window provided for distance to edge or midpoint calculation. ", "Add it using the addCell() function")
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
        out <- estPimsSingle(hypFrame[[x, "ppp"]], owins = hypFrame[x, "owins", drop = TRUE], pis = pis, null = null,
            tabObs = hypFrame[[x, "tabObs"]], centroids = hypFrame[x, "centroids", drop = TRUE], features = features,...)
        return(out)
    })
    list(hypFrame = hypFrame, null = null, pis = pis, features = features)
}
