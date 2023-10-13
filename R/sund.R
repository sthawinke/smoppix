#' Helper function to spit gene pairs
sund = function(x, sep = "_") {
    strsplit(x, sep)[[1]]
}
makeDesignVar = function(x, designVars, sep = "_"){
    apply(x[, designVars], 1, collapse = sep)
}
