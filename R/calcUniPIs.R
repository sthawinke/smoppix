#' Calculate univariate PIs for all features
#'
#' @inheritParams estPimsSingle
#' @param ecdfsCell ecdfs within the cells
#' @param ecdfAll Overall ecdf
#' @param nSub The number of events in the subsampled pattern
#' @param verbose A boolean, should verbose output be printed
#' @param cd A matrix of cross-distances
#'
#' @return PIs for every feature
#' @importFrom spatstat.geom marks
#' @importFrom BiocParallel bpparam
calcUniPIs <- function(p, pis, verbose, ecdfsCell, owins, tabObs, null, cd,
                       nPointsAll, nSub, ecdfAll, features, centroids) {
    if (verbose) {
        message("Calculating univariate probabilistic indices...")
    }
    splitFac <- rep(seq_len(bpparam()$workers), length.out = length(nams <- names(tabObs[features])))
    uniPIs <- unsplit(f = splitFac, bplapply(split(nams, f = splitFac), function(ss){
        lapply(ss, function(feat) {
            pSub <- p[id <- which(marks(p, drop = FALSE)$gene == feat), ]
            NP <- length(id)
            idLeft = if(noSubSam <- (nPointsAll >= (NPtot <- npoints(p)))){
                seq_len(NPtot)
            } else {
                sample(NPtot, nPointsAll)
            }
            obsDistNN2 = nndist(pSub)^2
            subCoords = getCoordsMat(pSub);
            leftCoords = getCoordsMat(if(noSubSam) p else p[idLeft])
            approxRanks = vapply(seq_along(id), FUN.VALUE = integer(1), function(i){
                findRanksDist(subCoords[i, ], leftCoords[idLeft != i,],
                              obsDistNN2[i])
            })
            ## TO DO: from here, calculate nnPI and alldistPI. No distance matrix needs to be stored
            # Avoid zero distances by removing observations of gene itself
            NNdistPI <- if (any(pis == "nn")) {
                if (NP == 1)
                  NA else calcNNPI(null, p = p, ecdfAll = ecdfAll, id = id)
            }
            # Also here room for improvement
            allDistPI <- if (any(pis == "allDist")) {
                if (NP == 1)
                  NA else calcAllDistPI(pSub, ecdfAll = ecdfAll, null = null, cd = cd)
            }
            rm(cd)
            edgeDistPI <- if (any(pis == "edge")) {
                calcWindowDistPI(pSub, owins, ecdfAll = ecdfsCell, pi = "edge")
            }
            midPointDistPI <- if (any(pis == "midpoint")) {
                calcWindowDistPI(pSub, owins, centroids, ecdfAll = ecdfsCell, pi = "midpoint")
            }
            nnCellPI <- if (any(pis == "nnCell")) {
                calcWindowDistPI(pSub, owins, ecdfAll = ecdfsCell, feat = feat, cellAndGene = marks(p, drop = FALSE)[,
                  c("gene", "cell")], pi = "nnCell", null = null, ecdfs = ecdfsCell)
            }
            allDistCellPI <- if (any(pis == "allDistCell")) {
                calcWindowDistPI(pSub, owins, ecdfAll = ecdfsCell, feat = feat, cellAndGene = marks(p, drop = FALSE)[,
                  c("gene", "cell")], pi = "allDistCell", null = null, ecdfs = ecdfsCell)
            }
            list(pointDists = c(nn = NNdistPI, allDist = allDistPI), windowDists = list(edge = edgeDistPI, nnCell = nnCellPI,
                allDistCell = allDistCellPI, midpoint = midPointDistPI))
        })
    }))
    names(uniPIs) <- nams
    return(uniPIs)
}
