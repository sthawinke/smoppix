#' A fast version that makes an ecdf based on an ordered vector and the cumulative probability function
#'
#' @param x The ordered vector
#' @param nuhe number of observations drawn from which the minuum is taken
#' @details
#'  A sampling trick: the negative hypergeometric distribution or the number of failures
#'  till the first success without replacement: the distribution of the minimum
#'
#' @return The ecdf function
ecdfFast = function (x, num)
{    rval <- approxfun(x, pnhyper(seq_along(x), m = num, n = length(x)-num, r = 1),
                      method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    rval
}
