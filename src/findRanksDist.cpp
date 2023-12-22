#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int findRanksDist(NumericVector m1, NumericMatrix m2, double squaredDist) {
    int nrow2 = m2.nrow();

    int out = 0;
    for (int r2 = 0; r2 < nrow2; r2++) {
       double dist = pow(m1(0) - m2(r2, 0), 2) + pow(m1(1) - m2(r2, 1), 2);
        if(dist < squaredDist){
            out ++;
        }
    }
    return out;
}
/*** R
a = matrix(rnorm(400), ncol = 2)
vec = rnorm(2)
dist = rnorm(1)
# Call new function
findRanksDist(vec, a, dist^2)
#sum(colSums((t(a)-vec)^2) < dist^2)
*/
