#' A fast version that makes an ecdf based on an presorted vector
#'
#' @param x The ordered vector
#' @return The ecdf function
#' @importFrom stats approxfun
ecdfPreSort = function (x)
{
    n <- length(x)
    rval <- approxfun(x, seq_along(x)/n,
                      method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    attr(rval, "call") <- sys.call()
    rval
}
