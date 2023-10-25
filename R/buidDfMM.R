#' Build a data frame for a single gene and measure to prepare for mixed model building
#'
#' @param pimRes The result of a call to estPims()
#' @param gene A character string indicating the desired gene
#' @param pi character string indicating the desired pim as outcome in the linear mixed model
#'
#' @return A dataframe
#' @export
#'
#' @examples
buildDfMM = function(pimRes, gene, pi = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint")){
    stopifnot(length(gene)==1)
    pi = match.arg(pi)
    piListName = if(pairId <- grepl("Pair", pi)) "biPIs" else "uniPIs"
    piListNameInner = if(winId <- any(pi == c("edge", "midpoint"))) "windowDists" else "pointDists"
    Design = names(pimRes)
    piEsts = lapply(pimRes, function(x){
        x[[piListName]][[gene]][[piListNameInner]][[pi]]
    })
}
