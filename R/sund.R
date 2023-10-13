#' Helper function to spit gene pairs
#'
#' @param x
#' @param sep
sund = function(x, sep = "_") {
    strsplit(x, sep)[[1]]
}
#' Make design variable by combining different design variables
#'
#' @param x
#' @param designVars
#' @param sep
#'
#' @return
makeDesignVar = function(x, designVars, sep = "_"){
    apply(x[, designVars], 1, collapse = sep)
}
