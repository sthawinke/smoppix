#' Extract all unique features from an object
#'
#' @param x A hyperframe or a results list containing a hyperframe
#'
#' @return A vector of features
getFeatures = function(x){
    x = if(is(x, "hyperframe")){
        x
    } else if(is(x$hypFrame, "hyperframe")){
        x$hypFrame
    }
    unique(unlist(lapply(x$tabObs, names)))
}
