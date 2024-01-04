#' Converta list of differently formatted windows to owins
#'
#' @param windows The list of windows
#'
#' @return A list of owins
convertToOwins = function(windows){
    out = if(all(vapply(windows, FUN.VALUE = TRUE, is.owin))){
        return(windows)
    }
}
