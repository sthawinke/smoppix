#' Subsample a point pattern when it is too large
#'
#' @param p The point pattern
#' @param nSims The maximum size
#'
#' @return A point pattern, subsample if necessary
subSampleP = function(p, nSims){
    if((NP <- npoints(p)) > nSims) p[sample(NP, nSims),] else p
}
