#' Calculate bivaraiet PIs
#'
#' @param manyPairs Number of gene pairs to be considered "many"
#' @param allowManyGenePairs A boolean, set to true to override warning on computation time
#' @inheritParams calcUniPIs
#'
#' @return A matrix of bivariate PIs
calcBiPIs = function(p, pis, null, nSims, ecdfs, nSub, ecdfAll, features,
                     manyPairs, verbose, allowManyGenePairs){
    if(!allowManyGenePairs &&
       (numGenePairs <- choose(length(features), 2)) > manyPairs){
        warning(immediate. = TRUE, "Calculating probablistic indices for",
                numGenePairs, "gene pairs may take a long time!\n",
                "Set allowManyGenePairs to TRUE to suppress this announcement")
    } else  if(verbose){
        message("Calculating bivariate probabilistic indices...")
    }
    genePairsMat = combn(features, 2)
    out = matrix(nrow = ncol(genePairsMat), byrow = TRUE,
                 vapply(seq_len(ncol(genePairsMat)), FUN.VALUE = double(length(grep(pis, pattern = "Pair"))),
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
