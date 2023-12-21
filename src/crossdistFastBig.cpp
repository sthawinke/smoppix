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
void crossdistFastInner(NumericMatrix m1, NumericMatrix m2, Mat<T> outBigMat) {
  int nrow1 = m1.nrow();
  int nrow2 = m2.nrow();
  int ncol = m1.ncol();

  if (ncol != m2.ncol()) {
    throw std::runtime_error("Incompatible number of dimensions");
  }

  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(m1(r1, c12) - m2(r2, c12), 2);
      }
      outBigMat(r1, r2) = sqrt(total);
    }
  }
}
// #Taken from https://www.r-bloggers.com/2019/09/matrix-cross-distances-using-rcpp/
// [[Rcpp::export]]
void crossdistFast(NumericMatrix m1, NumericMatrix m2, SEXP pOutBigMat) {
    // First we tell Rcpp that the object we've been given is an external
    // pointer.
    XPtr<BigMatrix> xpOutMat(pOutBigMat);
    crossdistFastInner(
        m1, m2,
        arma::Mat<double>((double *)xpOutMat->matrix(), xpOutMat->nrow(), xpOutMat->ncol(), false)
    );
    return;
    }

/*** R
a <- matrix(rnorm(16), ncol = 2)
b <- matrix(rnorm(18), ncol = 2)

# Call new euclidean function
bm_out <- spatrans:::crossdistWrapper(a, b, returnBigMatrix = TRUE)
*/
