#' Fit a gradient to a point pattern
#'
#' A Poisson process is fitted to the data assuming exponential relationship
#' between x and y variables and intensity
#'
#' @param hypFrame the hyperframe
#' @param returnModel A boolean, should the entire model be returned?
#' Otherwise the p-value is returned
#' @param ... passed onto spatstat.model::ppm
#' @inheritParams estGradients
#'
#' @return A list with p-value and coefficients,
#' or a mppm model when returnModel is true
#' @importFrom stats formula
#' @importFrom spatstat.model mppm anova.mppm
#' @seealso \link{estGradients}
fitGradient = function(hypFrame, fixedForm, randomForm, fixedFormSimple,
                       returnModel = FALSE, silent, ...){
    xyModel = try(mppm(data = hypFrame, fixedForm, random = randomForm, ...),
                  silent = silent)
     if(returnModel){
        return(xyModel)
    } else if(is(xyModel, "try-error")){
        return(list("pVal" = 1, "coef" = NULL))
   } else {
        intModel = try(mppm(data = hypFrame, fixedFormSimple, random = randomForm, ...), silent = silent)
        pVal = if(is(intModel, "try-error") ){
            NA
        } else {
            anovaRes = anova(xyModel, intModel, test = "Chisq")
            anovaRes$"Pr(>Chi)"[2]
        }
        out = list("pVal" = pVal, "coef" = coef(xyModel))
        return(out)
    }
}
