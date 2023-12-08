#' Check if the required pi's are present in the object
#'
#' @param x The result of the PI calculation, or a weight function
#' @param pi A character string indicating the desired PI
#'
#' @return Throws an error when the PIs are not found, otherwise returns invisible
checkAttr <- function(x, pi) {
  if (!(pi %in% attr(x, "pis"))) {
    stop(
      "Required PI not present in object. Rerun ",
      switch(class(x)[1],
        "list" = "estPims",
        "scam" = "buildWeightFunction"
      ),
      " with correct 'pis' argument"
    )
  } else {
    invisible()
  }
}
