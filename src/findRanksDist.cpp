#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix findRanksDist(NumericMatrix coords, NumericMatrix coordsLeft, NumericMatrix NNDist) {
    int nrow = coords.nrow();
    int nrowLeft = coordsLeft.nrow();
    int nDist = NNDist.ncol();

    IntegerMatrix out(nrow*2, nDist);
    //#First half for larger than, second for ties
    for (int r1 = 0; r1 < nrow; r1++) {//# Loop over events
        for (int r2 = 0; r2 < nrowLeft; r2++) {//# Loop over randomly drawn events
            double dist = sqrt(pow(coords(r1, 0) - coordsLeft(r2, 0), 2) +
                pow(coords(r1,1) - coordsLeft(r2, 1), 2));
            for (int r3 = 0; r3 < nDist; r3++) {//# Loop over distances
                if(dist < NNDist(r1, r3)){
                    out(r1, r3)++;
                } else if(dist == NNDist(r1, r3)){
                    out(r1+nrow, r3)++;
                    //#Keep track of ties too
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
out = findRanksDist(coordsMat, coordsMatLeft, distMat)
tmp = crossdist(X = coordsMat[, 1], Y = coordsMat[, 2], x2 = coordsMatLeft[, 1], y2 = coordsMatLeft[, 2])
out2 = vapply(seq_len(ncol(distMat)), FUN.VALUE = matrix(0,nrow(distMat), 2), function(x){
    cbind(rowSums(distMat[,x] > tmp), rowSums(distMat[,x] == tmp))
})
out2 = rbind(out2[,1,], out2[,2,])
all(out==out2)
*/
