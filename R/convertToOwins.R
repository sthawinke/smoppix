#' Convert windows in different formats to spatstat owin format
#' @description Convert a list of differently formatted windows to owins, for
#' addition to a hyperframe
#'
#' @param windows The list of windows
#' @param namePPP the name of the point pattern, will be added to the cell names
#' @inheritParams addCell
#' @details Order of traversion of polygons may differ between data types.
#' Where applicable, different orders are tried before throwing an error.
#'
#' @return A list of owins
#' @importFrom spatstat.geom is.owin owin as.owin
#' @importFrom methods slot is
convertToOwins <- function(windows, namePPP, coords, ...) {
    if (is(windows, "SpatialPolygonsDataFrame")) {
        if (requireNamespace("polyCub")) {
            p <- slot(windows, "polygons")
            winOut <- lapply(p, as.owin, ...)
        } else {
            stop("Install polyCub package first")
        }
    } else if (allClass(windows, is.owin)) {
        winOut <- windows
    } else if (allClass(windows, function(x) {
        is.data.frame(x) || is.matrix(x)
    })) {
        winOut <- lapply(windows, function(df) {
            i <- seq_len(nrow(df))
            foo <- try(silent = TRUE, owin(poly = list(
                x = df[, coords[1]],
                y = df[, coords[2]]
            ), ...))
            if (is(foo, "try-error")) {
                foo <- try(silent = TRUE, owin(poly = list(
                    x = df[rev(i), coords[1]],
                    y = df[, coords[2]]
                ), ...))
            }
            if (is(foo, "try-error")) {
                foo <- try(silent = TRUE, owin(poly = list(
                    x = df[, coords[1]],
                    y = df[rev(i), coords[2]]
                ), ...))
            }
            if (is(foo, "try-error")) {
                foo <- try(silent = TRUE, owin(poly = list(
                    x = df[rev(i), coords[1]],
                    y = df[rev(i), coords[2]]
                ), ...))
            }
            foo
        })
    } else if ((allRoi <- allClass(windows, is, "ijroi")) || is(windows, "ijzip")) {
        if (requireNamespace("RImageJROI")) {
            winOut <- if (allRoi) {
                lapply(windows, function(w){
                    foo <- try(silent = TRUE, RImageJROI::ij2spatstat(w))
                    if(is(foo, "try-error")){
                        w$coords = w$coords[rev(seq_len(nrow(w$coords))),]
                        foo <- try(silent = TRUE, RImageJROI::ij2spatstat(w))
                    }
                    return(foo)
                    }, ...)
            } else {
                RImageJROI::ij2spatstat(windows)
            }
        } else {
            stop("Install RImageJROI package first")
        }
    } else {
        stop("Input type unknown")
    }
    if (is.null(names(winOut))) {
        names(winOut) <- paste0("Cell", seq_along(winOut))
    }
    return(winOut)
}
allClass <- function(x, checkFun, ...) {
    all(vapply(x, FUN.VALUE = TRUE, checkFun, ...))
}
