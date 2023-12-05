#' Calculate the PI under the CSR null hypothesis,
#' weighing by Poisson variability
#'
#' @param NP The number of events
#' @param nAll The number of events in the calculation of the PI under the null
#' @param approxRanks The approximate ranks in the null distribution
#'
#' @return The PI
#' @importFrom stats qpois dpois
#' @importFrom extraDistr pnhyper
getPoissonPi = function(NP, nAll, approxRanks){
    poisQ = qpois(lambda = NP, p = c(1e-4, 1-1e-4))
    corFac = (1- sum(dpois(0:1, lambda = NP))) #Redistribute weights on 0 and 1
    sum(vapply(seq.int(max(poisQ[1], 2), poisQ[2]), FUN.VALUE = double(1), function(xx){
        mean(pnhyper(approxRanks, r = 1, m = xx - 1, n = nAll - xx + 1))*
            dpois(xx, lambda = NP)/corFac
    }))
}
