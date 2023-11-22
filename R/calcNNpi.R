#' Estimate the PI for the nearest neighbour distances
#'
#' @inheritParams estPimsSingle
#' @inheritParams calcAllDistPI
#'
#' @importFrom spatstat.geom nndist crossdist npoints area
#' @importFrom spatstat.random rpoispp
#' @importFrom stats ecdf
#' @importFrom Rfast rowMins colMins
#' @return The estimated probabilistic index
calcNNPI = function(pSub, p, null, nSims, rowSortMat){
    obsDistNN = nndist(pSub)
    if(null == "background"){
        NP = npoints(pSub)
        obsDistRank = rowMins(obsDistNN >= rowSortMat) #The ranks: first FALSE. Total number until first
        mean(pnhyper(obsDistRank, n = nrow(rowSortMat)-NP+1, m = NP-1, r = 1))#Exclude event itself
    } else if(null == "CSR"){
        lambda = npoints(pSub)/area(p$window)
        simDistsNN = lapply(integer(nSims), function(i){
            nndist(rpoispp(lambda, win = pSub$window))
        }) #Density dependent so still sims needed
        mean(ecdf(unlist(simDistsNN))(obsDistNN))
    }
}
calcNNPIpair = function(cd, id1, id2, null, p, rowSortMat){
    npp = npoints(p)
    if(null == "background"){
        obsDistNN = c(rowMins(cd, value = TRUE), colMins(cd, value = TRUE))
        obsDistRank = rowMins(obsDistNN >= rowSortMat[c(id1, id2),]) #The ranks: first FALSE
        seq1 = seq_along(id1);seq2 = seq_along(id2) + length(id1)
        mean(na.rm = TRUE, c(pnhyper(obsDistRank[seq1], n = npp-length(id2), m = length(id2), r = 1),
                             pnhyper(obsDistRank[seq2], n = npp-length(id1), m = length(id1), r = 1)))
        #Using the negative hypergeometric precludes Monte-Carlo
    } else if(null == "CSR"){
        pSim = runifpoint(npp, win = p$window)
        simDistsNN = lapply(integer(nSims), function(i){
            id <- sample(npp, length(id1)+length(id2), replace = Replace)
            p1 = pSim[id[seq_len(np1)],];p2 = pSim[-id[seq_len(np1)],]
            distMatSub = crossdist(p1, p2)
            c(rowMins(distMatSub, value = TRUE), colMins(distMatSub, value = TRUE))
        })
        mean(ecdf(unlist(simDistsNN))(obsDistNN))
    }
}
#' Fast version of nncross
#'
#' @param p1,p2 The point patterns
#'
#' @return A vector of nearest neighbour distances
#' @importFrom Rfast rowMins
#' @importFrom spatstat.geom crossdist
nncrossFast = function(p1, p2){
    rowMins(crossdist(p1,p2), value = TRUE)
}
