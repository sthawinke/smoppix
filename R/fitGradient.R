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
#' @return A vector of of length 3 with p-value and variances of the random effects,
#' or a mppm model when returnModel is true
#' @importFrom stats formula
#' @importFrom spatstat.model mppm anova.mppm
#' @importFrom nlme VarCorr
#' @seealso \link{estGradients}
fitGradient = function(hypFrame, fixedEffects = NULL, randomEffects = NULL, returnModel = FALSE, silent,...){
    fixedForm = buildFormula(outcome = "ppp", fixedVars = c(fixedEffects, "x:id", "y:id", "id"),
                             randomVars = NULL)
    randomForm = buildFormula(outcome = "", fixedVars = NULL, randomVars = randomEffects)
    xyModel = try(mppm(data = hypFrame, fixedForm, random = randomForm, verb = FALSE, ...),
                  silent = silent)
     if(returnModel){
        return(xyModel)
    } else if(is(xyModel, "try-error")){
        return(c("pVal" = 1, "x" = NA, "y" = NA))
   } else {
       fixedFormReduced = buildFormula(outcome = "ppp", fixedVars = c(fixedEffects, "id"),
                                randomVars = NULL)
        intModel = try(mppm(data = hypFrame, fixedFormReduced, random = randomForm, verb = FALSE, ...), silent = silent)
        pVal = if(is(intModel, "try-error") ){
            NA
        } else {
            anovaRes = anova(xyModel, intModel, test = "Chisq")
            if(anovaRes["xyModel", "logLik"] < anovaRes["intModel", "logLik"]){
                1 #If elaborate model not better, return 1
            } else {anovaRes$"Pr(>Chi)"[2]}
        }
        out = list("pVal" = pVal, "coef" = coef(xyModel))
        return(out)
    }
}
