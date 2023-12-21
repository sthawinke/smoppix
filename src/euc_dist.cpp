// To enable the functionality provided by Armadillo's various macros,
// simply include them before you include the RcppArmadillo headers.
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

// The following header file provides the definitions for the BigMatrix
// object
#include <bigmemory/BigMatrix.h>

// C++11 plugin
// [[Rcpp::plugins(cpp11)]]

template <typename T>
void BigArmaEuclidean(const Mat<T>& inBigMat, Mat<T> outBigMat) {

    int W = inBigMat.n_rows;

    for(int i = 0; i < W - 1; i++){
        for(int j=i+1; j < W; j++){
            outBigMat(j-1,i) = sqrt(sum(pow((inBigMat.row(i) - inBigMat.row(j)),2)));
        }
    }
}

// [[Rcpp::export]]
void BigArmaEuc(SEXP pInBigMat, SEXP pOutBigMat) {
    // First we tell Rcpp that the object we've been given is an external
    // pointer.
    XPtr<BigMatrix> xpMat(pInBigMat);
    XPtr<BigMatrix> xpOutMat(pOutBigMat);
    int type = xpMat->matrix_type();
    switch(type) {
    case 1:
        BigArmaEuclidean(
            arma::Mat<char>((char *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
            arma::Mat<char>((char *)xpOutMat->matrix(), xpOutMat->nrow(), xpOutMat->ncol(), false)
        );
        return;

    case 2:
        BigArmaEuclidean(
            arma::Mat<short>((short *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
            arma::Mat<short>((short *)xpOutMat->matrix(), xpOutMat->nrow(), xpOutMat->ncol(), false)
        );
        return;

    case 4:
        BigArmaEuclidean(
            arma::Mat<int>((int *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
            arma::Mat<int>((int *)xpOutMat->matrix(), xpOutMat->nrow(), xpOutMat->ncol(), false)
        );
        return;

    case 8:
        BigArmaEuclidean(
            arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
            arma::Mat<double>((double *)xpOutMat->matrix(), xpOutMat->nrow(), xpOutMat->ncol(), false)
        );
        return;

    default:
        // We should never get here, but it resolves compiler warnings.
        throw Rcpp::exception("Undefined type for provided big.matrix");
    }

}


/*** R
library(bigmemory)
mat <- matrix(rnorm(16), 4)
bm <- as.big.matrix(mat)

# Call new euclidean function
bm_out <- spatrans:::bigMatrixEuc(bm)
*/
