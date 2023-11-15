#' Helper function to spit gene pairs
#'
#' @param x character string
#' @param sep The character used to split
sund = function(x, sep = "--") {
    strsplit(x, sep)[[1]]
}
#' Make design variable by combining different design variables
#'
#' @param x the design matrix
#' @param designVars the design variables to be combined
#' @param sep The string to separate the components
#'
#' @return a vector of design levels
makeDesignVar = function(x, designVars, sep = "_"){
    apply(x[, designVars], 1, collapse = sep)
}
