#' Fit a gradient to a point pattern
#'
#' A Poisson process is fitted to the data assuming exponential relationship
#' between x and y variables and intensity
#'
#' @param ppp The point pattern
#' @param trend A formula, describing the log intensity
#' @param returnModel A boolean, should the entire model be returned?
#' Otherwise the gradient is returned
#' @param ... passed onto spatstat.model::ppm
#'
#' @return A ppm object, or a scalar that is the gradient
#' @importFrom stats formula
#' @importFrom spatstat.model ppm coef.ppm
#' @seealso \link{estGradients}
fitGradient = function(ppp, trend = formula("~ x + y"), returnModel = FALSE, ...){
    mod = ppm(ppp, trend = trend, ...)
    if(returnModel){
        return(mod)
    } else {
        Coef = coef.ppm(mod)[, c("x", "y")]
        Grad = sqrt(Coef["x"]^2 + Coef["y"]^2)
        return(Grad)
    }
}
