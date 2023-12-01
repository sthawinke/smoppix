#' Remove the diagonal from a square matrix
#' @details The diagonal elements are removed, and the elements right from the
#' diagonal are shifted left to form a new matrix with one column less
#'
#' @param x
#'
#' @return A matrix with the same number of rows as x and one column less
dropDiagonal = function(x){
    diagId = seq(1, length(x), by = (NR <- nrow(x))+1)
    t(matrix(t(x)[-diagId], ncol = NR))
}
