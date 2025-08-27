#' Parallel processing with BiocParallel with load balancing
#' @description The vector to iterate over (iterator) is split into as many parts as there are
#' cores available, such that each core gets an equal load and overhead is minimized. 
#' The registered backend is then used by default to multithread using \link[BiocParallel]{bplapply}.
#'
#' @param iterator The vector to iterate over
#' @param func The function to apply to each element
#' @param loopFun The looping function, can also be 'lapply' for serial processing
#'
#' @return A list with the same length as iterator
#' @importFrom BiocParallel bpparam bplapply bpnworkers
#' @export
loadBalanceBplapply <- function(iterator, func, loopFun = "bplapply") {
    loopFunMatched <- match.fun(loopFun)
    if (loopFun == "lapply") {
        loopFunMatched(iterator, func)
    } else if (loopFun == "bplapply") {
        splitFac <- rep(seq_len(bpnworkers(bpparam())), length.out = length(iterator))
        # Divide the work over the available workers
        unsplit(f = splitFac, loopFunMatched(
            split(iterator, f = splitFac),
            function(ss) {
                out <- lapply(ss, func)
                names(out) <- ss
                return(out)
            }
        ))
    }
}
