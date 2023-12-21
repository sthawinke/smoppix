#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix crossdistFastLocal(NumericMatrix m1, NumericMatrix m2) {
    int nrow1 = m1.nrow();
    int nrow2 = m2.nrow();
    int ncol = m1.ncol();

    if (ncol != m2.ncol()) {
        throw std::runtime_error("Incompatible number of dimensions");
    }

    NumericMatrix out(nrow1, nrow2);

    for (int r1 = 0; r1 < nrow1; r1++) {
        for (int r2 = 0; r2 < nrow2; r2++) {
            double total = 0;
            for (int c12 = 0; c12 < ncol; c12++) {
                total += pow(m1(r1, c12) - m2(r2, c12), 2);
            }
            out(r1, r2) = sqrt(total);
        }
    }

    return out;
}
// #Taken from https://www.r-bloggers.com/2019/09/matrix-cross-distances-using-rcpp/

