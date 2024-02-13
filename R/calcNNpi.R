#' Estimate the PI for the nearest neighbour distances, given a set of ranks,
#' using the negative hypergeometric distribution
#'
#' @param n the total number of observed distances minus the number of distances
#' under consideration (the number of black balls in the urn)
#' @param Ranks The (approximate) ranks
#' @param m the number of observed distances (white balls in the urn)
#' @param r The rank of distances considered, r=1 is nearest neighbour distance
#'
#' @importFrom extraDistr pnhyper
#' @return A vector of evaluations of the negative hypergeometric distribution
#' function
#' @seealso \link{pnhyper}, \link{calcIndividualPIs}
calcNNPI <- function(Ranks, n, m, r = 1) {
    pnhyper(Ranks, n = n, m = m, r = r)
    # Exclude event itself, so NP - 1 m = N-1: White balls, number of other
    # events of the same gene n = n - (NP-1): black balls, number of events of
    # other genes in background r=1: Nearest neighbour so first occurrence
}
