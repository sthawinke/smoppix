#' Fit a gradient to a point pattern
#'
#' A Poisson process is fitted to the data assuming exponential relationship
#' between x and y variables and intensity
#'
#' @param hypFrame the hyperframe
#' @param fixedEffects Character vector of fixed effects present in the hyperframe,
#' modifying the baseline intensity. See details.
#' @param returnModel A boolean, should the entire model be returned?
#' Otherwise the p-value is returned
#' @param ... passed onto spatstat.model::ppm
#' @details The fixed effects modify the baseline intensity of the point pattern, not the gradient!
#'
#' @return A vector of of length 3 with p-value and variances of the random effects,
#' or a mppm model when returnModel is true
#' @importFrom stats formula
#' @importFrom spatstat.model mppm anova.mppm
#' @importFrom nlme VarCorr
#' @seealso \link{estGradients}
fitGradient = function(hypFrame, fixedEffects = NULL, returnModel = FALSE, silent,...){
    fixedForm = paste("ppp ~ 1", paste(paste("+",  fixedEffects), collapse = ""))
    xyModel = try(mppm(data = hypFrame, fixedForm, random = "x + y |id", ...), silent = silent)
     if(returnModel || is(xyModel, "try-error")){
        return(xyModel)
    } else {
        intModel = try(mppm(data = hypFrame, fixedForm, random = "1 |id", ...), silent = TRUE)
        pVal = if(is(intModel, "try-error")){
            NA
        } else {
            anovaRes = anova(xyModel, intModel, test = "Chisq")
            anovaRes$p.value[2]
        }
        Variances = as.numeric(VarCorr(xyModel$Fit$FIT)[c("x", "y"), "Variance"])
        out = c(pVal, Variances);names(out) = c("pVal", "x", "y")
        return(out)
    }
}
