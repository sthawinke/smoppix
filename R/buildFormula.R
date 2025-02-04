#' Build a formula from different components
#'
#' @param Formula A formula. If not supplied or equals NULL, will be overridden
#' @param fixedVars,randomVars Character vectors with fixed and random variables
#' @param outcome A character vector describing the outcome
#'
#' @return A formula
#' @details Random intercepts are assumed for the random effects, if more
#' complicated designs are used, do supply your own formula.
#' @importFrom stats formula
#' @seealso \link{fitLMMs},\link[stats]{formula}
buildFormula <- function(Formula, fixedVars, randomVars, outcome = "pi - 0.5") {
    # Allow cell as design variable, both fixed and random
    Formula <- if (missing(Formula) || is.null(Formula)) {
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
characterFormula <- function(Formula) {
    paste(as.character(Formula)[c(2, 1, 3)], collapse = " ")
}
