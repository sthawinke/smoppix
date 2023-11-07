#' Check if the required pi's are present in the object
#'
#' @param pimRes The result of the PI calculation
#' @param pi A character string indicating the desired PI
#'
#' @return Throws an error when the PIs are not found, otherwise returns invisible
checkAttr = function(pimRes, pi){
    if(!(pi %in% attr(pimRes, "pis"))){
        stop("Required PI not present in object. Rerun estPims with correct 'pis' argument")
    } else {invisible()}
}
