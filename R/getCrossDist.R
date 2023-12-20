#' Get a cross distance matrix, avoiding zero entries
#'
#' @param p The total point pattern
#' @param pSub A smaller sampled point pattern, representing the null
#' @param id The elements of p present in pSub
#' @inheritParams estPims
#'
#' @return A cross distance matrix without zeroes
#' @importFrom spatstat.geom npoints
getCrossDist <- function(p, pSub, null, id = NULL) {
    cd <- crossDistProxy(p, pSub)
    nSub <- npoints(pSub)
    if (null == "background") {
        if (npoints(p) == nSub) {
            # If all points used, drop diagonal
            cd <- dropDiagonal(cd)
        } else {
            if(is.null(id)){
                idDupl <- which(idMat <- (cd == 0), arr.ind = TRUE)
              # Points resampled, and thus distance 0. Replace by other points
                cd[idMat] <- crossDistProxy(p[idDupl[, 1], ],
                            p[sample(seq_len(npoints(p))[-idDupl[, 1]], 1), ])
            } else {
                diag(cd[id,]) = crossDistProxy(p[id, ],
                                    p[sample(seq_len(npoints(p))[-id], 1), ])
            }
        }
    }
    return(cd)
}
#' @importFrom proxy dist
#' @importFrom spatstat.geom coords
crossDistProxy = function(p, pSub){
    proxy::dist(coords(p), coords(pSub))
}
