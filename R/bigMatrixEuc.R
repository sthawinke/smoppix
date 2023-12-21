bigMatrixEuc <- function(bigMat){
    zeros <- big.matrix(nrow = nrow(bigMat)-1,
                        ncol = nrow(bigMat)-1,
                        init = 0,
                        type = typeof(bigMat))
    BigArmaEuc(bigMat@address, zeros@address)
    return(zeros)
}
