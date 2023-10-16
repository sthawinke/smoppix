#' Check in which cell each event lies and add a cell marker
#'
#' @param hypFrame The hyperframe
#' @param rois the list of list of rois per hyperframe
#' @return the modified hyperframe
addCell = function(hypFrame, rois){
    stopifnot(length(hypFrame) == length(rois), all(names(hypFrame) %in% names(rois)))
}
