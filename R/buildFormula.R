#' Build a formula from different components
#'
#' @param Formula A formula. If not supplied (is null), will be overridden
#' @param fixedVars,randomVars Character vectors with fixed and random variables
#' @param outcome A character vector describing the outcome
#'
#' @return A formula
#' @details Random intercepts are assumed for the random effects, if more
#' complicated designs are used, do supply your own formula
#' @importFrom stats formula
buildFormula = function(Formula, fixedVars, randomVars, outcome = "pi - 0.5"){
    # Allow cell as design variable, both fixed and random
    Formula =  if (is.null(Formula)) {
        fixedPart <- paste(outcome, " ~ 1", if (!is.null(fixedVars)) {
            paste("+", paste(fixedVars, collapse = " + "))
        })
        Formula <- formula(formChar <- paste(fixedPart, if (length(randomVars)) {
            paste("+", paste0("(1|", randomVars, ")", collapse = " + "))
        }, collapse = " "))
    } else {
        formula(Formula)
    }
    return(Formula)
}
characterFormula = function(Formula){
    paste(as.character(Formula)[c(2,1,3)], collapse = " ")
}
#' @importFrom stats formula
getFixedPart = function(Formula){
    if(!any(grepl("\\|", Formula))){
        Formula
    } else {
    formula(paste(collapse = " ", grep(value = TRUE, invert = TRUE, "\\|",
                        strsplit(as.character(Formula), split = "\\+")[[1]])))
    }
}
