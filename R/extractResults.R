#' Extract results from a list of fitted LMMs
#'
#' @param models The models
#' @param fixedVars The fixed effects for which the effect is to be reported
#' @param method passed onto p.adjust
#' @importFrom stats p.adjust
#' @import lmerTest
#' @return A list of matrices, all containing estimate, standard error, p-value and ajdusted p-value
#' @seealso \link{fitLMMs}
extractResults = function(models, fixedVars = NULL, method = "BH"){
    ints = t(vapply(models, FUN.VALUE = double(3), function(x){
        if(is.null(x) || is(x, "try-error")){
            c("Estimate" = NA, "Std. Error" = NA, "Pr(>|t|)" = NA)
        } else {
        summary(x)$coef["(Intercept)", c("Estimate", "Std. Error", "Pr(>|t|)")]
        }
    }))
    colnames(ints) = c("Estimate", "SE", "pVal")
    ints[, "Estimate"] = ints[, "Estimate"] + 0.5
    intMat = cbind(ints, "pAdj" = p.adjust(ints[, "pVal"], method = method))[order(ints[, "pVal"]),]
    #Order by p-value
    AnovaTabs = lapply(models[id <- vapply(models, FUN.VALUE = TRUE, is, "lmerModLmerTest")], anova)
    fixedOut = lapply(fixedVars, function(Var){
        pVal = vapply(AnovaTabs, FUN.VALUE = double(1), function(x) x[Var, "Pr(>F)"])
        coefs = lapply(models[id], function(model){
            coefObj = summary(model)$coef
            Coefs = coefObj[grep(paste0(Var, "[[:digit:]]"), rownames(coefObj)), "Estimate"]
        })
        coefMat = matrix(unlist(coefs), byrow = TRUE, nrow = length(pVal))
        colnames(coefMat) = paste0(Var, seq_len(ncol(coefMat)))
        cbind(coefMat, "pVal" = pVal,
              "pAdj" = p.adjust(pVal, method = method))[order(pVal),]
    })
    names(fixedOut) = fixedVars
    list("Intercept" = intMat, "fixedEffects" = fixedOut)
}

