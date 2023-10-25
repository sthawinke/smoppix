#' Try to find overlap between list of windows, and throw error when found
#'
#' @param owins the list of windows
#'
#' @return Throws an error when overlap found, otherwise returns invisible
#' @importFrom utils combn
findOverlap = function(owins){
    combs = combn(length(owins), 2)
    lapply(seq_len(ncol(combs)), function(i){
        if(overlap.owin(owins[[combs[1, i]]], owins[[combs[2, i]]])>0){
            stop("Overlap detected between windows ", combs[1, i], " and ", combs[2, i],
                 "!\nWindows must be non-overlapping.")
        }
    })
    return(invisible())
}
