#' Converta list of differently formatted windows to owins
#'
#' @param windows The list of windows
#'
#' @return
#' @export
#'
#' @examples
convertToOwins = function(windows){
    if(all(vapply(windows, FUN.VALUE = TRUE, is.owin))){
        return(windows)
    }
}
