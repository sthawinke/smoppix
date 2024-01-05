#' Converta list of differently formatted windows to owins
#'
#' @param windows The list of windows
#' @inheritParams addCell
#'
#' @return A list of owins
#' @importFrom spatstat.geom is.owin owin
#' @importFrom methods slot is
convertToOwins = function(windows, coords){
    if(allClass(windows, is.owin)){
        return(windows)
    } else if(allClass(windows, is, "SpatialPolygonsDataFrame")){
        if(require(polyCub)){
            lapply(windows, function(spdf){
                p <- slot(spdf, "polygons")
                lapply(p, polyCub::as.owin.SpatialPolygons)
            })
        } else {
            stop("Install polyCub package first")
        }
    } else if(allClass(windows, function(x) {is.data.frame(x) || is.matrix(x)})){
        lapply(windows, function(df){
            i = seq_len(nrow(df))
            foo = try(silent = TRUE, owin(poly = list(x = df[, coords[1]],
                                                      y = df[, coords[2]])))
            if(is(foo, "try-error")){
                foo = try(silent = TRUE, owin(poly = list(x = df[rev(i), coords[1]],
                                                          y = df[, coords[2]])))
            }
            if(is(foo, "try-error")){
                foo = try(silent = TRUE, owin(poly = list(x = df[, coords[1]],
                                                          y = df[rev(i), coords[2]])))
            }
            if(is(foo, "try-error")){
                foo = try(silent = TRUE, owin(poly = list(x = df[rev(i), coords[1]],
                                                          y = df[rev(i), coords[2]])))
            }
            foo
        })
    }
}
allClass  = function(x, checkFun, ...){
    all(vapply(x, FUN.VALUE = TRUE, checkFun, ...))
}
