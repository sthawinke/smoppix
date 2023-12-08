#' A wrapper function for the different probabilistic indices (PIs)
#'
#' @param p The point pattern
#' @param pis The probabilistic indices to be estimated
#' @param null A character vector, indicating how the null distribution is defined. See details.
#' @param nPointsAll How many points to subsample or simulate to calculate overall interpoint distance
#' and distance to point
#' @param nPointsAllWin How many points to subsample or simulate to calculate distance to cell edge or midpoint distribution
#' @param allowManyGenePairs A boolean, set to true to suppress warning messages for large numbers of gene pairs
#' @param manyPairs An integer, what are considered many gene pairs
#' @param verbose Should verbose output be printed?
#' @param owins,centroids The list of windows corresponding to cells, and their centroids
#' @param features A character vector, for which features should the probabilistic indices be calculated?
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
estPimsSingle <- function(p, pis, null, tabObs, nPointsAll = 5e2,
                          nPointsAllWin = 2e2, features = NULL,
                          allowManyGenePairs = FALSE, manyPairs = 1e6, verbose = FALSE,
                          owins = NULL, centroids = NULL) {
    if (is.null(features)) {
        features <- names(tabObs)
    } else {
        features <- intersect(features, names(tabObs))
    }
    if (null == "CSR") {
        if (any(pis %in% c("allDist", "allDistPair", "nn", "nnPair"))) {
            pSim <- runifpoint(nPointsAll, win = p$window)
            nSub <- npoints(pSim)
        }
        if (any(pis %in% c("allDist", "allDistPair", "nn", "nnPair"))) {
            ecdfAll <- ecdf(dist(coords(pSim)))
        }
        if (any(pis %in% c(
            "edge", "midpoint", "nnCell", "allDistCell",
            "nnPairCell", "allDistPairCell"
        ))) {
            ecdfsCell <- lapply(names(owins), function(nam) {
                pSub <- runifpoint(nPointsAllWin, win = owins[[nam]])
                edge <- if (any(pis == "edge")) {
                    ecdf(nncross(pSub, edges(owins[[nam]]), what = "dist"))
                }
                midpoint <- if (any(pis == "midpoint")) {
                    ecdf(crossdist(pSub, centroids[[nam]]))
                }
                allDistCell <- if (any(pis %in% c(
                    "nnCell", "allDistCell",
                    "nnPairCell", "allDistPairCell"
                ))) {
                    ecdf(dist(coords(pSub)))
                }
                list(
                    "edge" = edge, "midpoint" = midpoint,
                    "allDistCell" = allDistCell
                )
            })
            names(ecdfsCell) <- names(owins)
        }
    } else if (null == "background") {
        if (any(pis %in% c("edge", "midpoint"))) {
            ecdfsCell <- lapply(names(owins), function(nam) {
                pSub <- subSampleP(
                    pCell <- p[marks(p, drop = FALSE)$cell == nam, ],
                    nPointsAllWin
                )
                edge <- if (any(pis == "edge")) {
                    ecdf(nncross(pSub, edges(owins[[nam]]), what = "dist"))
                }
                midpoint <- if (any(pis == "midpoint")) {
                    ecdf(crossdist(pSub, centroids[[nam]]))
                }
                allDistCell <- if (any(pis %in% c(
                    "nnCell", "allDistCell",
                    "nnPairCell", "allDistPairCell"
                ))) {
                    apply(
                        rowSort(getCrossDist(pCell, pSub, null)), 1,
                        ecdfPreSort
                    )
                }
                list(
                    "edge" = edge, "midpoint" = midpoint,
                    "allDistCell" = allDistCell
                )
            })
            names(ecdfsCell) <- names(owins)
        }
        if (any(pis %in% c("nn", "nnPair", "allDist", "allDistPair"))) {
            # Prepare some null distances, either with CSR or background
            subSamSize <- max(c(nPointsAll, tabObs + 1e2))
            pSub <- subSampleP(p, subSamSize)
            # Subsample for memory reasons
            cd <- getCrossDist(p, pSub, null)
            ecdfs <- apply(rowSort(cd), 1, ecdfPreSort)
            # For background, condition on point locations
            nSub <- ncol(cd)
        }
    }
    # Univariate patterns
    if (any(idOne <- (tabObs[features] == 1)) && verbose) {
        message(
            "Features\n", paste(collapse = ", ", names(tabObs[features])[idOne]),
            "\nhave only one observation, so no intervent distances were calculated for them\n"
        )
    }
    piPair <- grepl(pis, pattern = "Pair")
    uniPIs <- if (any(!piPair)) {
        calcUniPIs(p, pis, verbose, ecdfsCell, owins, tabObs, null,
            ecdfs, nSub, ecdfAll, idOne, features,
            centroids = centroids
        )
    }
    # Bivariate patterns
    biPIs <- if (any(piPair)) {
        calcBiPIs(
            features = features[tabObs[features] > 0], p, pis, null,
            ecdfs, nSub, ecdfAll, manyPairs, verbose, allowManyGenePairs,
            ecdfsCell
        )
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
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPims(hypYang, pis = c("nn", "nnPair"))
#' # Both univariate and bivariate tests, with nearest neighbour distances
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
#' The suffix "Cell" indicates distances are being calculated within cells only.
estPims <- function(hypFrame, pis = c(
                        "nn", "allDist", "nnPair", "allDistPair",
                        "edge", "midpoint", "nnCell", "allDistCell",
                        "nnPairCell", "allDistPairCell"
                    ),
                    null = c("background", "CSR"),
                    features = attr(hypFrame, "features"), ...) {
    pis <- match.arg(pis, several.ok = TRUE)
    null <- match.arg(null)
    Im <- attr(hypFrame, "imageVars")
    if (any(pis %in% c("edge", "midpoint", "nnCell", "allDistCell")) && is.null(hypFrame$owins)) {
        stop(
            "No window provided for distance to edge or midpoint calculation.",
            "Add it using the addCell() function"
        )
    }
    if (any(id <- !(features %in% attr(hypFrame, "features")))) {
        stop("Features ", features[id], " not found in hyperframe")
    }
    hypFrame$pimRes <- bplapply(seq_len(nrow(hypFrame)), function(x) {
        with(hypFrame[x, , drop = TRUE], estPimsSingle(ppp,
            owins = owins,
            pis = pis, null = null, tabObs = tabObs, centroids = centroids, ...
        ))
    })
    attr(hypFrame, "pis") <- pis # Tag the pims calculated
    attr(hypFrame, "featuresEst") <- features
    attr(hypFrame, "imageVars") <- Im
    # Remember for which features the pims were calculated
    hypFrame
}
