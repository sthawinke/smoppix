library(testthat)
library(spatrans)
library(spatstat.random)
library(BiocParallel)
n <- 1e3 # number of molecules
ng <- 25 # number of genes
nfov = 3 # Number of fields of view
conditions = 3
# sample xy-coordinates in [0, 1]
x <- runif(n)
y <- runif(n)
# assign each molecule to some gene-cell pair
gs <- paste0("gene", seq(ng))
gene <- sample(gs, n, TRUE)
fov <- sample(nfov, n, TRUE)
condition = sample(conditions, n, TRUE)
# construct data.frame of molecule coordinates
df <- data.frame(gene, x, y, fov, "condition" = condition)
#A list of point patterns
listPPP = tapply(seq(nrow(df)), df$fov, function(i){
    ppp(x = df$x[i], y = df$y[i], marks = df[i, c("gene", "condition"), drop = FALSE])
}, simplify = FALSE)
# Regions of interest (roi): Diamond in the center plus four triangles
w1 <- owin(poly=list(x=c(0,.5,1,.5),y=c(.5,0,.5,1)))
w2 <- owin(poly=list(x=c(0,0,.5),y = c(.5,0,0)))
w3 <- owin(poly=list(x=c(0,0,.5),y=c(1,0.5,1)))
w4 <- owin(poly=list(x=c(1,1,.5),y=c(0.5,1,1)))
w5 <- owin(poly=list(x=c(1,1,.5),y=c(0,0.5,0)))
wWrong <- owin(poly=list(x=c(0,1,1,0),y=c(.75,.5,.75,1)))
hypFrame <- buildHyperFrame(df, coordVars = c("x", "y"), imageVars = c("condition", "fov"))
nDesignFactors = length(unique(hypFrame$image))
wList = lapply(seq_len(nDesignFactors), function(x){
    list("w1" = w1, "w2" = w2, "w3" = w3, "w4" = w4, "w5" = w5)
})
wList2 = lapply(seq_len(nDesignFactors), function(x){
    list("w1" = w1, "w2" = w2, "w3" = w3, "w4" = w4, "w5" = w5, "wWrong" = wWrong)
})
names(wList) = names(wList2) = rownames(hypFrame)
hypFrame2 = addCell(hypFrame, wList)
#Register the parallel backend
#register(SerialParam())
nCores = 2
register(MulticoreParam(nCores))
piEstsBG <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint", "midpoint"),
                    point = c(0.5, 0.5), null = "background")
piEstsCSR <- estPims(hypFrame2, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "fixedpoint", "midpoint"),
                     point = c(0.5, 0.5), null = "CSR")
objBG <- addWeightFunction(piEstsBG, designVars = "condition")
objCSR <- addWeightFunction(piEstsCSR, designVars = "condition")
test_check("spatrans")
