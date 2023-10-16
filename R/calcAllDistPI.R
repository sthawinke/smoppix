calcAllDistPI = function(pSub, ecdfAll){
    obsDist = pairdist(pSub)
    mean(ecdfAll(obsDist))
}
