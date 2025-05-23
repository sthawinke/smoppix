#' Estimate the PI for the nearest neighbour distances, given a set of ranks,
#' using the negative hypergeometric distribution
#'
#' @param Ranks The (approximate) ranks, number of times observed distance is larger
#' @param n the total number of observed distances minus the number of distances
#' under consideration (the number of failures or black balls in the urn)
#' @param m the number of observed distances (successes or white balls in the urn)
#' @param r The rank of distances considered, r=1 is nearest neighbour distance
#' @param ties The number of times the observed distance is equal to a null distance,
#' of the same length as Ranks
#' @details Ties are counted half to match the definition of the PI.
#' @importFrom extraDistr pnhyper
#' @return A vector of evaluations of the negative hypergeometric distribution
#' function
#' @seealso \link[extraDistr]{pnhyper}, \link{calcIndividualPIs}
calcNNPI <- function(Ranks, n, m, ties, r = 1) {
    tmp <- pnhyper(Ranks, n = n, m = m, r = r)
    if (!is.null(ties)) {
        id <- ties > 0
        tmp[id] <- 0.5 * (tmp[id] + pnhyper(Ranks[id] + ties[id], n = n, m = m, r = r))
        # Everything up to Ranks counted twice, up to ties+Ranks counted once,
        # so divide by 2
    }
    # Exclude event itself, so NP - 1 m = N-1: White balls, number of other
    # events of the same gene n = n - (NP-1): black balls, number of events of
    # other genes in background r=1: Nearest neighbour so first occurrence
    # Equalities are counted half
    return(tmp)
}
