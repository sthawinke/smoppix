# Helper function to get gene pair or order reversed
getGp = function(x, gp){
    if(isVec <- (is.vector(x) || is.list(x))){
        Names = names(x)
    } else if(isMat <- is.matrix(x)){
        Names = colnames(x)
    }
    if(gp %in% Names) {
        if(isVec) x[[gp]] else x[, gp]
    } else {
        geneSplit = strsplit(gp, "_")[[1]]
        gp = paste(rev(geneSplit), collapse = "_")
        if(gp %in% Names) {
            if(isVec) x[[gp]] else x[, gp]
        } else {NULL}
    }
}
