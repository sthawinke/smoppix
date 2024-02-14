#' Fit a gradient to a point pattern
#'
#' A Poisson process is fitted to the data assuming exponential relationship
#' between x and y variables and intensity
#'
#' @param ppp The point pattern
#' @param trend A formula, describing the log intensity
#' @param returnModel A boolean, should the entire model be returned?
#' Otherwise the gradient and its variance are returned
#' @param ... passed onto spatstat.model::ppm
#'
#' @return A ppm object, or a scalar that is the gradient
#' @importFrom stats formula
#' @importFrom spatstat.model ppm coef.ppm vcov.ppm
#' @importFrom spatstat.geom unmark is.marked
#' @seealso \link{estGradients}
fitGradient = function(ppp, trend = formula("~ x + y"), returnModel = FALSE, ...){
    if(is.marked(ppp)){
        ppp = unmark(ppp)
        #Get rid of marks for ppm.ppp function
    }
    mod = try(ppm(ppp, trend = trend, ...))
    if(returnModel){
        return(mod)
    } else if(is(mod, "try-error")){
        return(c("gradient" = NA, "variance" = NA))
    } else {
        Coef = coef.ppm(mod)[c("x", "y")]
        Grad = sqrt(Coef["x"]^2 + Coef["y"]^2) #The gradient
        Deriv = Coef/Grad #Derivative to x- and y-slopes
        VarGrad = if(is.null(vcov.ppm(mod))){
            NA
        } else {
           crossprod(Deriv, vcov.ppm(mod)[c("x", "y"), c("x", "y")]) %*% Deriv
        }
        #When Fisher information matrix is singular, no variance can be found,
        #so attach no weight
        return(c("gradient" = unname(Grad), "variance" = c(VarGrad)))
    }
}
