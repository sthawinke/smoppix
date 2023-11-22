#' An aux function to build gene pairs
#'
#' @param genes The genes to be combined
#'
#' @return A chcracter vector of gene pairs
makePairs = function(genes){
    apply(combn(genes, 2), 2, paste, collapse = "--")
}
