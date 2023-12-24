#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix findRanksDist(NumericMatrix coords, NumericMatrix coordsLeft, NumericMatrix squaredNNDist) {
    int nrow = coords.nrow();
    int nrowLeft = coordsLeft.nrow();
    int nDist = squaredNNDist.ncol();

    IntegerMatrix out(nrow, nDist);
    for (int r1 = 0; r1 < nrow; r1++) {//# Loop over events
        for (int r2 = 0; r2 < nrowLeft; r2++) {//# Loop over randomly drawn events
            double dist = pow(coords(r1, 0) - coordsLeft(r2, 0), 2) +
                pow(coords(r1,1) - coordsLeft(r2, 1), 2);
            for (int r3 = 0; r3 < nDist; r3++) {//# Loop over distances
                if(squaredNNDist(r1, r3) > dist){
                    out(r1, r3)++;
                }
            }
        }
    }
    return out;
}
/*** R
coordsMat = matrix(rnorm(400), ncol = 2)
coordsMatLeft = matrix(rnorm(600), ncol = 2)
distMat = matrix(rnorm(1200), 200, 6)
# Call new function
out = findRanksDist(coordsMat, coordsMatLeft, distMat^2)
library(proxy)
tmp = proxy::dist(coordsMat, coordsMatLeft)^2
out2 = vapply(seq_len(ncol(distMat)), FUN.VALUE = double(nrow(distMat)), function(x){
    rowSums(distMat[,x]^2 > tmp)
})
*/
