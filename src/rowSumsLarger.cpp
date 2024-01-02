#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// The following header file provides the definitions for the BigMatrix
// object
#include <bigmemory/BigMatrix.h>
// C++11 plugin
// [[Rcpp::plugins(cpp11)]]

template <typename T>
IntegerVector rowSumsL(NumericVector m1, Mat<T> inBigMat) {
    const int len1 = m1.length();
    const int ncol = inBigMat.n_cols;

    IntegerVector out(len1);

    for (int r1 = 0; r1 < len1; r1++) {
        int total = 0;
        for (int r2 = 0; r2 < ncol; r2++) {
                if(m1(r1) > inBigMat(r1, r2)){
                    total += 1;
                }
        }
        out(r1) = total;
    }
    return out;
}

// [[Rcpp::export]]
IntegerVector rowSumsLarger(NumericVector m1, SEXP pInBigMat) {
    // First we tell Rcpp that the object we've been given is an external
    // pointer.
    XPtr<BigMatrix> xpInMat(pInBigMat);
    IntegerVector out = rowSumsL(
        m1, arma::Mat<double>((double *)xpInMat->matrix(), xpInMat->nrow(), xpInMat->ncol(), false)
    );
    return out;
}

/*** R
library(bigmemory)
a <- as.big.matrix(matrix(rnorm(32), ncol = 8, nrow = 4))
b <- rnorm(4)
# Call new function
bm_out <- spatrans:::rowSumsWrapper(a, b)
#bm_out_check <- spatrans:::rowSumsWrapper(a[], b)
*/