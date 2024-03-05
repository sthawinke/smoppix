#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix findRanksDist(NumericMatrix coords, NumericMatrix coordsLeft, NumericMatrix squaredNNDist) {
    int nrow = coords.nrow();
    int nrowLeft = coordsLeft.nrow();
    int nDist = squaredNNDist.ncol();

    NumericMatrix out(nrow, nDist);
    for (int r1 = 0; r1 < nrow; r1++) {//# Loop over events
        for (int r2 = 0; r2 < nrowLeft; r2++) {//# Loop over randomly drawn events
            double dist = pow(coords(r1, 0) - coordsLeft(r2, 0), 2) +
                pow(coords(r1,1) - coordsLeft(r2, 1), 2);
            for (int r3 = 0; r3 < nDist; r3++) {//# Loop over distances
                if(squaredNNDist(r1, r3) > dist){
                    out(r1, r3)++;
                } else if(squaredNNDist(r1, r3) == dist){
                    out(r1, r3) = out(r1, r3) + 0.5;
                }
            }
        }
    }
    return out;
}
/*** R
coordsMat = matrix(rpois(400, 2), ncol = 2)
coordsMatLeft = matrix(rpois(600, 3), ncol = 2)
distMat = matrix(rpois(1200, 4), 200, 6)
# Call new function
out = findRanksDist(coordsMat, coordsMatLeft, distMat^2)
tmp = crossdist(X = coordsMat[, 1], Y = coordsMat[, 2], x2 = coordsMatLeft[, 1], y2 = coordsMatLeft[, 2])^2
out2 = vapply(seq_len(ncol(distMat)), FUN.VALUE = double(nrow(distMat)), function(x){
    rowSums(distMat[,x]^2 > tmp) + 0.5*rowSums(distMat[,x]^2 == tmp)
})
all(out==out2)
*/
