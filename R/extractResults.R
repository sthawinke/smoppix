#' Extract results from a list of fitted LMMs
#'
#' @param models The models
#' @param fixedVars The fixed effects for which the effect is to be reported
#' @param method passed onto p.adjust
#' @importFrom stats p.adjust
#' @return A list of matrices, all containing estimate, standard error, p-value and ajdusted p-value
#' @seealso \link{fitLMMs}
extractResults = function(models, fixedVars, method = "BH"){
    out = lapply(Vars <- c("(Intercept)", fixedVars), function(Var){
        if(Var == "(Intercept)"){
            ints = vapply(models, FUN.VALUE = double(3), function(x){
                summary(x)$coef["(Intercept)", c("Estimate", "Std. Error", "Pr(>|t|)")]
               })
            colnames(tmp) = c("Estimate", "SE", "pVal")
            tmp[, "Estimate"] = tmp[, "Estimate"] + 0.5
            cbind(tmp, "pAdj" = p.adjust(tmp[, "pVal"], method = method))[order(tmp[, "pVal"]),]
            #Order by p-value
        } else {
            Coefs = lapply(models, function(x){
                summary(x)$coef[grep(Var, rownames(summary(x)$coef)),
                                c("Estimate", "Std. Error", "Pr(>|t|)")]
            })
        }
    })
    names(out) = Vars
    out
}
