#' Calculate bivaraiet PIs
#'
#' @param manyPairs Number of gene pairs to be considered 'many'
#' @param allowManyGenePairs A boolean, set to true to override warning on computation time
#' @param maxNum The maximum number of points to calculate cross distances for,
#' to avoid excessive memory usage.
#' @inheritParams calcUniPIs
#' @importFrom spatstat.geom nncross
#' @importFrom BiocParallel bpparam
#'
#' @return A matrix of bivariate PIs
calcBiPIs <- function(p, pis, null, cd, nSub, ecdfAll, features, manyPairs, verbose, allowManyGenePairs, ecdfsCell, nPointsAll,
    maxNum = 10000) {
    if (!allowManyGenePairs && (numGenePairs <- choose(length(features), 2)) > manyPairs) {
        warning(immediate. = TRUE, "Calculating probablistic indices for", numGenePairs, "gene pairs may take a long time!\n",
            "Set allowManyGenePairs to TRUE to suppress this announcement")
    } else if (verbose) {
        message("Calculating bivariate probabilistic indices...")
    }
    genePairsMat <- combn(features, 2)
    splitFac <- rep(seq_len(bpparam()$workers), length.out = ncol(genePairsMat))
    tmp <- unsplit(f = splitFac, bplapply(split(seq_len(ncol(genePairsMat)), f = splitFac), function(ss) {
        lapply(ss, function(i) {
            feat1 <- genePairsMat[1, i]
            feat2 <- genePairsMat[2, i]
            pSub1 <- p[id1 <- which(marks(p, drop = FALSE)$gene == feat1), ]
            pSub2 <- p[id2 <- which(marks(p, drop = FALSE)$gene == feat2), ]
            cd = try(silent = TRUE, crossdistWrapper(p[c(id1, id2),], subSampleP(p[-c(id1, id2),], nPointsAll),
                 returnBigMatrix = prod(length(id1) + length(id2), nPointsAll) > 1e4))
            if(is(cd, "try-error")){
                gc()
                cd = crossdistWrapper(p[c(id1, id2),], subSampleP(p[-c(id1, id2),], nPointsAll),
                                  returnBigMatrix = prod(length(id1) + length(id2), nPointsAll) > 1e4)
            }
            # Reorder and subset if needed
            NNdistPI <- if (any(pis == "nnPair")) {
                calcNNPIpair(obsDistNN = c(nncross(pSub1, pSub2, what = "dist"),
                                           nncross(pSub2, pSub1, what = "dist")),
                             id1 = id1, id2 = id2, null = null,
                             cd = cd, ecdfAll = ecdfAll, n = nSub)#[c(id1, id2),]
            }
            nnCellPI <- if (any(pis == "nnPairCell")) {
                calcWindowPairPI(pSub1, pSub2, ecdfAll = ecdfsCell, ecdfs = ecdfsCell, pi = "nnPairCell", null = null,
                  feat1 = feat1, feat2 = feat2, cellAndGene = marks(p, drop = FALSE)[, c("gene", "cell")])
            }
            list(pointDists = c(nnPair = NNdistPI), windowDists = list(nnPairCell = nnCellPI))
        })
    }))
    names(tmp) <- apply(genePairsMat, 2, paste, collapse = "--")
    return(tmp)
}
