#' Estimate probabilistic indices for single-molecule localization patterns, and add the variance weighting function.
#' @description Estimate different probabilistic indices for localization
#' on all point patterns of a hyperframe, and integrate the results in the same hyperframe. 
#' estPisSingle() is the workhorse function for a single point pattern.
#' @export
#' @param hypFrame A hyperframe
#' @param ... additional arguments, passed on to \link{estPisSingle}.
#' @param verbose A boolean, whether to report on progress of the fitting process.
#' @param null A character vector, indicating how the null distribution is
#'  defined. See details.
#' @param nPointsAll,nPointsAllWithinCell How many points to subsample or simulate to calculate the
#' overall nearest neighbour distance distribution under the null hypothesis.
#' The second argument (nPointsAllWithinCell) applies to within cell calculations, where a lower number usually suffises.
#' @param nPointsAllWin How many points to subsample or simulate
#' to calculate distance to cell edge or centroid distribution
#' @param minDiff An integer, the minimum number of events from other genes
#'  needed for calculation of background distribution of distances.
#'  Matters mainly for within-cell calculations: cells with too few events are skipped.
#' @param minObsNN An integer, the minimum number of events required for a gene to be analysed. See details.
#' @param features A character vector, for which features should the
#' probabilistic indices be calculated?
#' @return For estPis(), the hyperframe with the estimated PIs present in it
#' @order 1
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' yangPims <- estPis(hypYang[c(seq_len(4), seq(27, 29)), ], pis = "nn",
#'     nPointsAll = 4e2)
#' # Univariate nearest neighbour distances
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root")) 
#' # Add the weight functions 
#' yangObj <- addWeightFunction(yangPims, lowestLevelVar = "section", 
#' pi = "nn")
#' # Alternative formulation with 'lowestLevelVar'
#' @details
#' The null distribution used to calculate the PIs can be either 'background' or 'null'.
#' For 'background', the observed distributions of all genes is used. Alternatively,
#' for null = 'CSR', Monte-Carlo simulation under complete spatial randomness
#' is performed within the given window to find the null distribution of the distance under study.
#' See \insertCite{Hawinkel2025b}{smoppix} for precise definition of the PI.
#'
#' The 'nn' prefix indicates that nearest neighbour distances are being used,
#' either univariately or bivariately.
#' The suffix 'Pair' indicates that bivariate probabilistic indices,
#' testing for co- and antilocalization are being used.
#' 'edge' and 'centroid' calculate the distance to the edge respectively
#'  the centroid of the windows added using the \link{addCell} function.
#' The suffix 'Cell' indicates that nearest neighbour distances are being calculated within cells only.
#'
#' It can be useful to set the minObsNN higher than the default of 5 for calculations within cells when the number of events is low,
#' not to waste computation time on gene (pairs) with very variable PI estimates.
estPis <- function(
    hypFrame, pis = c("nn", "nnPair", "edge", "centroid", "nnCell",
        "nnPairCell"), verbose = TRUE, null = c("background", "CSR"),
    nPointsAll = switch(null,
        background = 5e4, CSR = 2e3
    ), nPointsAllWithinCell = switch(null,
        background = 5e3, CSR = 1e3
    ), nPointsAllWin = 1e4, minDiff = 20, minObsNN = 1L,
    features = getFeatures(hypFrame), ...) {
    pis <- match.arg(pis, several.ok = TRUE)
    null <- match.arg(null)
    stopifnot(
        is.hyperframe(hypFrame), is.numeric(nPointsAll), is.numeric(nPointsAllWithinCell),
        is.numeric(nPointsAllWin), is.numeric(minDiff), is.character(features)
    )
    features <- sort(features)
    if (any(vapply(hypFrame$ppp, FUN.VALUE = FALSE, function(x) is.unsorted(getCoordsMat(x)[, "x"])))) {
        stop("Point patterns must be sorted by x-coordinates!
        Build the hyperframe by the buildHyperFrame function or presort!")
    }
    if (any(pis %in% c("edge", "centroid", "nnCell")) && is.null(hypFrame$owins)) {
        stop(
            "No window provided for distance to edge or centroid calculation. ",
            "Add it using the addCell() function"
        )
    }
    foo <- checkFeatures(hypFrame, features)
    if (verbose) {
        message("Calculating PIs for point pattern ")
    }
    hypFrame$pimRes <- lapply(seq_len(nrow(hypFrame)), function(x) {
        if (verbose) {
            message(x, " of ", nrow(hypFrame), "  ")
        }
        out <- estPisSingle(hypFrame[[x, "ppp"]],
            owins = hypFrame[x, "owins", drop = TRUE],
            pis = pis, null = null, tabObs = hypFrame[[x, "tabObs"]], centroids = hypFrame[x,
                "centroids",
                drop = TRUE
            ], features = features, nPointsAll = nPointsAll,
            nPointsAllWithinCell = nPointsAllWithinCell, nPointsAllWin = nPointsAllWin,
            minDiff = minDiff, minObsNN = minObsNN, ...
        )
        return(out)
    })
    list(hypFrame = hypFrame, null = null, pis = pis, features = features)
}
#' @param p The point pattern
#' @param owins,centroids The list of windows corresponding to cells,
#' and their centroids
#' @param tabObs A table of observed gene frequencies
#' @param window An window of class owin, in which events can occur
#' @param loopFun The function to use to loop over the features.
#' Defaults to bplapply except when looping over features within cells
#'
#' @return For estPisSingle(), a list of data frames with estimated PIs per gene and/or gene pair:
#' \item{pointDists}{PIs for pointwise distances overall}
#' \item{windowDists}{PIs for distances to cell wall or centroid}
#' \item{withinCellDists}{PIs for pointwise distances within cell}
#' @importFrom stats ecdf dist
#' @importFrom spatstat.random runifpoint
#' @importFrom spatstat.geom is.hyperframe
#' @importFrom Rdpack reprompt
#' @rdname estPis
#' @order 2
estPisSingle <- function(
    p, pis, null, tabObs, owins = NULL, centroids = NULL, window = p$window,
    loopFun = "bplapply", features, nPointsAll, nPointsAllWithinCell, nPointsAllWin,
    minDiff, minObsNN) {
    features <- intersect(features, names(tabObs))
    if(!length(features)){
        return(list(pointDists = NULL, windowDists = NULL, withinCellDists = NULL))
    }
    # Scramble to ensure equal calculation times in multithreading
    if (any(pis %in% c("nn", "nnPair"))) {
        if (null == "background") {
            # Subset some background points
            pSubLeft <- subSampleP(p, nPointsAll, returnId = TRUE)
        } else if (null == "CSR") {
            pSim <- runifpoint(nPointsAll, win = window)
            ecdfAll <- ecdf(stats::dist(getCoordsMat(pSim)))
        }
    }
    if (any(pis %in% c("edge", "centroid"))) {
        ecdfsCell <- findEcdfsCell(
            p = p, owins = owins, centroids = centroids, nPointsAllWin = nPointsAllWin,
            null = null, pis = pis, loopFun = loopFun
        )
    }
    piList <- if (!all(idCell <- grepl("Cell", pis))) {
        calcIndividualPIs(
            p = p, pSubLeft = pSubLeft, pis = pis, null = null, tabObs = tabObs,
            owins = owins, centroids = centroids, features = features, ecdfAll = ecdfAll,
            ecdfsCell = ecdfsCell, minDiff = minDiff, loopFun = loopFun, minObsNN = minObsNN
        )
    }
    nnPis <- if ("nn" %in% pis) {
        vapply(piList[intersect(features, names(tabObs[tabObs > 1]))],
            FUN.VALUE = double(1),
            function(x) {
                x$pointDists$nn
            }
        )
    }
    if ("nnPair" %in% pis && length(features) > 1) {
        genePairsMat <- combn(features, 2)
        nnPairPis <- vapply(seq_len(ncol(genePairsMat)), FUN.VALUE = double(1), function(i) {
            feat1 <- genePairsMat[1, i]
            feat2 <- genePairsMat[2, i]
            mean(c(getElement(piList[[feat1]]$pointDists$nnPair, feat2), 
                   getElement(piList[[feat2]]$pointDists$nnPair,feat1))
                 )
            #Each point gets the same weight, variance weighting is implemented downstream
        })
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
                p = pSub, null = null, nPointsAll = nPointsAllWithinCell, window = owins[[nam]],
                features = features, tabObs = table(marks(pSub, drop = FALSE)$gene),
                loopFun = "lapply", minDiff = minDiff, minObsNN = minObsNN
            )$pointDists
        }, loopFun = loopFun)
        names(cellDists) <- unCells
        cellDists
    }
    list(pointDists = pointDists, windowDists = windowDists, withinCellDists = withinCellDists)
}
