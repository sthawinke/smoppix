#' Calculate the Moran's I test statistic for spatial autocorrelation
#'
#' The normalized Moran's I test statistic is calculated, with standard normal
#' distribution under the null hypothesis
#' @note The implementation is inspired on the one from ape::Moran.I, but more
#'  bare-bones faster and using less memory.
#'
#' @param x A vector of outcomes
#' @param W The matrix of weights, with dimensions equal to the lenght of x
#'
#' @return The normalized Moran's I statistic
moransI = function(x, W){
    n = length(x)
    stopifnot(n==ncol(W))
    rawCen = x-mean(x) #Mean centered raw data
    v = sum(rawCen^2)
    mI = (tcrossprod(rawCen, W) %*% rawCen)/v + 1/(n-1)
    #Moran's I, centered around zero by 1/(n-1), the expected value
    S1 <- 0.5 * sum((W + t(W))^2)
    S2 <- sum((1 + colSums(W))^2)
    s.sq <- n^2
    k <- (sum(rawCen^4)/n)/(v/n)^2
    sdi <- sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq) -
                     k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq))/((n -1) *
                                                                          (n - 2) * (n - 3) * s.sq) - 1/((n - 1)^2))
    return(as.double(mI/sdi))
}
