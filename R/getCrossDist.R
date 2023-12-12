#' Get a cross distance matrix, avoiding zero entries
#'
#' @param p The total point pattern
#' @param pSub A smaller sampled point pattern, representing the null
#' @inheritParams estPims
#'
#' @return A cross distance matrix without zeroes
#' @importFrom spatstat.geom npoints crossdist
getCrossDist <- function(p, pSub, null) {
  cd <- crossdist(p, pSub)
  nSub <- npoints(pSub)
  if (null == "background") {
    if (npoints(p) == nSub) {
      # If all points used, drop diagonal
      cd <- dropDiagonal(cd)
    } else {
      idDupl <- which(idMat <- (cd == 0), arr.ind = TRUE)
      # Points resampled, and thus distance 0. Replace by other points
      cd[idMat] <- crossdist(
        p[idDupl[, 1], ],
        p[sample(seq_len(npoints(p))[-idDupl], 1), ]
      )
    }
  }
  return(cd)
}
