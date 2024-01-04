#' Converta list of differently formatted windows to owins
#'
#' @param windows The list of windows
#' @inheritParams addCell
#'
#' @return A list of owins
#' @importFrom spatstat.geom is.owin owin
#' @importFrom polyCub as.owin.SpatialPolygons
#' @importFrom methods slot is
convertToOwins = function(windows, coords){
    if(allClass(windows, is.owin)){
        return(windows)
    } else if(allClass(windows, is, "SpatialPolygonsDataFrame")){
        if(require(polyCub)){
            lapply(windows, function(spdf){
                p <- slot(spdf, "polygons")
                lapply(p, as.owin.SpatialPolygons)
            })
        } else {
            stop("Install polyCub package first")
        }
    } else if(allClass(windows, function(x) {is.data.frame(x) || is.matrix(x)})){
        lapply(windows, function(df){
            foo = try(silent = TRUE, owin(poly = list(x = roisText[i, "X"],
                                                      y = roisText[i, "Y"])))
            if(is(foo, "try-error")){
                foo = try(silent = TRUE, owin(poly = list(x = roisText[rev(i), "X"],
                                                          y = roisText[i, "Y"])))
            }
            if(is(foo, "try-error")){
                foo = try(silent = TRUE, owin(poly = list(x = roisText[i, "X"],
                                                          y = roisText[rev(i), "Y"])))
            }
            if(is(foo, "try-error")){
                foo = try(silent = TRUE, owin(poly = list(x = roisText[rev(i), "X"],
                                                          y = roisText[rev(i), "Y"])))
            }
            foo
        })
    }
}
allClass  = function(x, checkFun, ...){
    all(vapply(x, FUN.VALUE = TRUE, checkFun, ...))
}
