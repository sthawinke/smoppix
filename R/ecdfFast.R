#' A fast version that makes an ecdf based on an ordered vector and the cumulative probability function
#'
#' @param x The ordered vector
#' @param cumDens Corresponding cumulative density
#'
#' @return The ecdf function
ecdfFast = function (x, cumDens)
{
    rval <- approxfun(x, cumDens,
                      method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    rval
}
