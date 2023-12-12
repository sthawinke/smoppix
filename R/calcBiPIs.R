#' Calculate bivaraiet PIs
#'
#' @param manyPairs Number of gene pairs to be considered "many"
#' @param allowManyGenePairs A boolean, set to true to override warning on computation time
#' @inheritParams calcUniPIs
#'
#' @return A matrix of bivariate PIs
calcBiPIs <- function(p, pis, null, ecdfs, nSub, ecdfAll, features,
                      manyPairs, verbose, allowManyGenePairs, ecdfsCell) {
  if (!allowManyGenePairs &&
    (numGenePairs <- choose(length(features), 2)) > manyPairs) {
    warning(
      immediate. = TRUE, "Calculating probablistic indices for",
      numGenePairs, "gene pairs may take a long time!\n",
      "Set allowManyGenePairs to TRUE to suppress this announcement"
    )
  } else if (verbose) {
    message("Calculating bivariate probabilistic indices...")
  }
  genePairsMat <- combn(features, 2)
  tmp <- lapply(seq_len(ncol(genePairsMat)), function(i) {
    feat1 <- genePairsMat[1, i]
    feat2 <- genePairsMat[2, i]
    pSub1 <- p[id1 <- which(marks(p, drop = FALSE)$gene == feat1), ]
    pSub2 <- p[id2 <- which(marks(p, drop = FALSE)$gene == feat2), ]
    cd <- crossdist(pSub1, pSub2)
    # Reorder and subset if needed
    NNdistPI <- if (any(pis == "nnPair")) {
      calcNNPIpair(
        cd = cd, id1 = id1, id2 = id2, null = null,
        ecdfs = ecdfs[c(id1, id2)], p = p,
        ecdfAll = ecdfAll, n = nSub
      )
    }
    allDistPI <- if (any(pis == "allDistPair")) {
      calcAllDistPIpair(
        id1 = id1, id2 = id2, ecdfAll = ecdfAll,
        null = null, ecdfs = ecdfs[c(id1, id2)],
        crossDist = cd
      )
    }
    nnCellPI <- if (any(pis == "nnPairCell")) {
      calcWindowPairPI(pSub1, pSub2,
        ecdfAll = ecdfsCell,
        ecdfs = ecdfsCell, pi = "nnPairCell", null = null,
        feat1 = feat1, feat2 = feat2, cd = cd,
        cellAndGene = marks(p, drop = FALSE)[, c("gene", "cell")]
      )
    }
    allDistCellPI <- if (any(pis == "allDistPairCell")) {
      calcWindowPairPI(pSub1, pSub2,
        cd = cd, null = null,
        feat1 = feat1, feat2 = feat2,
        cellAndGene = marks(p, drop = FALSE)[, c("gene", "cell")],
        ecdfAll = ecdfsCell, pi = "allDistPairCell",
        ecdfs = ecdfsCell
      )
    }
    list(
      "pointDists" = c("nnPair" = NNdistPI, "allDistPair" = allDistPI),
      "windowDists" = list("allDistPairCell" = allDistCellPI,
                           "nnPairCell" = nnCellPI)
    )
  })
  names(tmp) <- apply(genePairsMat, 2, paste, collapse = "--")
  return(tmp)
}
