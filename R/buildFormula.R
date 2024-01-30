buildFormula = function(Formula, fixedVars, randomVars, outcome = "pi - 0.5"){
    # Allow cell as design variable, both fixed and random
    if (is.null(Formula)) {
        fixedPart <- paste(outcome, " ~ 1", if (!is.null(fixedVars)) {
            paste("+", paste(fixedVars, collapse = " + "))
        })
        Formula <- formula(formChar <- paste(fixedPart, if (!is.null(randomVars)) {
            paste("+", paste0("(1|", randomVars, ")", collapse = " + "))
        }))
    }
    return(Formula)
}
characterFormula = function(Formula){
    paste(as.character(Formula), collapse = "")
}
