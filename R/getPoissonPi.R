#' Calculate the PI under the CSR null hypothesis,
#' weighing by Poisson variability
#'
#' @param NP The number of events
#' @param nAll The number of events in the calculation of the PI under the null
#' @param approxRanks The approximate ranks in the null distribution
#' @param pair A boolean, is a gene pair involved
#' @param quants The quantile of the Poisson distribution calculated by
#'
#' @return The PI
#' @importFrom stats qpois dpois
#' @importFrom extraDistr pnhyper
getPoissonPi = function(NP, nAll, approxRanks, pair = FALSE, quants = c(1e-2, 1-1e-2)){
    poisQ = qpois(lambda = NP, p = quants)
    Max = if(pair) 1 else 2  #For a gene pair, go down to 1 and do not exclude point itself
    corFac = (1- sum(dpois(0:Max, lambda = NP))) #Redistribute weights on 0 and 1
    sum(vapply(seq.int(max(poisQ[1], Max), poisQ[2]), FUN.VALUE = double(1), function(xx){
        mean(pnhyper(approxRanks, r = 1, m = xx -(!pair), n = nAll - xx + (!pair)))*
            dpois(xx, lambda = NP)/corFac
    }))

}
