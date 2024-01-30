#' A wrapper function for the different probabilistic indices (PIs)
#'
#' @param p The point pattern
#' @param owins,centroids The list of windows corresponding to cells,
#' and their centroids
#' @param tabObs A table of observed gene frequencies
#' @param window An window of class owin, in which events can occur
#' @param loopFun The function to use to loop over the features.
#' Defaults to bplapply except when looping over features within cells
#' @inheritParams estPis
#'
#' @return A list of data frames with estimated PIs per gene and/or gene pair:
#' \item{pointDists}{PIs for pointwise distances overall}
#' \item{windowDists}{PIs for distances to cell wall or centroid}
#' \item{withinCellDists}{PIs for pointwise distances within cell}
#' @importFrom stats ecdf dist
#' @importFrom spatstat.random runifpoint
#' @importFrom Rdpack reprompt
estPisSingle <- function(p, pis, null, tabObs, owins = NULL, centroids = NULL, window = p$window, loopFun = "bplapply",
                          features, nPointsAll, nPointsAllWithinCell, nPointsAllWin, minDiff ) {
    if (!length(features)) {
        return(list(pointDists = NULL, windowDists = NULL, withinCellDists = NULL))
    }
    features <- sample(intersect(features, names(tabObs)))
    # Scramble to ensure equal calculation times in multithreading
    if (any(pis %in% c("nn", "nnPair"))) {
        if (null == "background") {
            # Subset some background points
            pSubLeft <- subSampleP(p, nPointsAll, returnId = TRUE)
        } else if (null == "CSR") {
            pSim <- runifpoint(nPointsAll, win = window)
            ecdfAll <- ecdf(stats::dist(coords(pSim)))
        }
    }
    if (any(pis %in% c("edge", "centroid"))) {
        ecdfsCell <- findEcdfsCell(
            p = p, owins = owins, centroids = centroids,
            nPointsAllWin = nPointsAllWin, null = null, pis = pis,
            loopFun = loopFun
        )
    }
    piList <- if(!all(idCell <- grepl("Cell", pis))){
        calcIndividualPIs(
        p = p, pSubLeft = pSubLeft, pis = pis,
        null = null, tabObs = tabObs, owins = owins,
        centroids = centroids, features = features, ecdfAll = ecdfAll,
        ecdfsCell = ecdfsCell, minDiff = minDiff,
        loopFun = loopFun
    )}
    nnPis <- if ("nn" %in% pis) {
        vapply(piList[intersect(features, names(tabObs[tabObs > 1]))],
            FUN.VALUE = double(1), function(x) {
                x$pointDists$nn
            }
        )
    }
    if ("nnPair" %in% pis && length(features) > 1) {
        genePairsMat <- combn(features, 2)
        nnPairPis <- vapply(seq_len(ncol(genePairsMat)),
            FUN.VALUE = double(1),
            function(i) {
                feat1 <- genePairsMat[1, i]
                feat2 <- genePairsMat[2, i]
                mean(c(
                    getElement(piList[[feat1]]$pointDists$nnPair, feat2),
                    getElement(
                        piList[[feat2]]$pointDists$nnPair,
                        feat1
                    )
                ))
            }
        )
        names(nnPairPis) <- apply(genePairsMat, 2, paste, collapse = "--")
    } else {
        nnPairPis <- NULL
    }
    pointDists <- list(nn = nnPis, nnPair = nnPairPis)
    # Window pis
    windowDists <- if (any(pis %in% c("edge", "centroid"))) {
        lapply(piList, function(x) x$windowDists)
    }
    # Within cell: recurse into estPisSingle but now per cell
    withinCellDists <- if (any(idCell)) {
        unCells <- setdiff(unique(marks(p, drop = FALSE)$cell), "NA")
        cellDists <- loadBalanceBplapply(unCells, function(nam) {
            pSub <- p[marks(p, drop = FALSE)$cell == nam, ]
            estPisSingle(
                pis = gsub("Cell", "", grep(value = TRUE, "Cell", pis)),
                p = pSub, null = null,
                nPointsAll = nPointsAllWithinCell, window = owins[[nam]], features = features,
                tabObs = table(marks(pSub, drop = FALSE)$gene),
                loopFun = "lapply", minDiff = minDiff
            )$pointDists
        })
        names(cellDists) <- unCells
        cellDists
    }
    list(
        pointDists = pointDists, windowDists = windowDists,
        withinCellDists = withinCellDists
    )
}
#' Estimate probabilistic indices for (co)localization on a hyperframe
#' @description Estimate different probabilistic indices for (co)localization
#' on  a hyperframe, and integrate the results in the same hyperframe
#' @export
#' @param hypFrame the hyperframe
#' @param ... additional arguments, passed on to estPisSingle
#' @param verbose a boolean, whether to report on progress of fitting
#' @param pis The probabilistic indices to be estimated
#' @param null A character vector, indicating how the null distribution is
#'  defined. See details.
#' @param nPointsAll,nPointsAllWithinCell How many points to subsample or simulate to calculate
#' overall interpoint distance and distance to point.
#' This parameter is a strong driver of computation time.
#' The second argument applies to within cell calculations.
#' @param nPointsAllWin How many points to subsample or simulate
#' to calculate distance to cell edge or centroid distribution
#' @param minDiff An integer, the minimum number of events from other genes
#'  needed for calculation of background PIs
#' @param features A character vector, for which features should the
#' probabilistic indices be calculated?
#' @return The hyperframe with the estimated PIMs present in it
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPis(hypYang[c(seq_len(5), seq(25, 29)), ], pis = "nn")
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
#' 'edge' and 'centroid' calculate the distance to the edge respectively
#'  the centroid of the windows added using the addCell() function.
#' The suffix 'Cell' indicates distances are being calculated within cells only.
estPis <- function(hypFrame, pis = c("nn", "nnPair", "edge", "centroid", "nnCell", "nnPairCell"),
                    verbose = TRUE, null = c("background", "CSR"),
                    nPointsAll = switch(null, background = 1e4, CSR = 1e3),
                    nPointsAllWithinCell = switch(null, background = 1e4, CSR = 5e2
                    ), nPointsAllWin = 1e3,
                    minDiff = 2e1, features = getFeatures(hypFrame), ...) {
    pis <- match.arg(pis, several.ok = TRUE)
    null <- match.arg(null)
    if (any(pis %in% c("edge", "centroid", "nnCell")) && is.null(hypFrame$owins)) {
        stop("No window provided for distance to edge or centroid calculation. ", "Add it using the addCell() function")
    }
    if (any(id <- !(features %in% getFeatures(hypFrame)))) {
        stop("Features ", features[id], " not found in hyperframe")
    }
    if (verbose) {
        message("Calculating PIs for point pattern ")
    }
    hypFrame$pimRes <- lapply(seq_len(nrow(hypFrame)), function(x) {
        if (verbose) {
            message(x, " of ", nrow(hypFrame), "  ")
        }
        out <- estPisSingle(hypFrame[[x, "ppp"]],
            owins = hypFrame[x, "owins", drop = TRUE], pis = pis, null = null,
            tabObs = hypFrame[[x, "tabObs"]], centroids = hypFrame[x, "centroids", drop = TRUE],
            features = features, nPointsAll = nPointsAll, nPointsAllWithinCell = nPointsAllWithinCell,
            nPointsAllWin = nPointsAllWin, minDiff = minDiff, ...
        )
        return(out)
    })
    list(hypFrame = hypFrame, null = null, pis = pis, features = features)
}