library(testthat)
library(spatrans)
library(spatstat.random)
library(BiocParallel)
n <- 1e3 # number of molecules
ng <- 50 # number of genes
cellsnc <- 20 # number of
nfov = 3 # Number of fields of view
# sample xy-coordinates in [0, 1]
x <- runif(n)
y <- runif(n)
# assign each molecule to some gene-cell pair
gs <- paste0("gene", seq(ng))
gene <- sample(gs, n, TRUE)
fov <- sample(nfov, n, TRUE)
# assure gene & cell are factors so that
# missing observations aren't dropped
gene <- factor(gene, gs)
# construct data.frame of molecule coordinates
df <- data.frame(gene, x, y, fov)
#A list of point patterns
listPPP = tapply(seq(nrow(df)), df$fov, function(i){
    ppp(x = df$x[i], y = df$y[i], marks = df[i, "gene", drop = FALSE])
}, simplify = FALSE)
# Regions of interest (roi)
w1 <- owin(poly=list(x=c(0,.25,.5,.25),y=c(.75,.5,.75,1)))
w2 <- owin(poly=list(x=c(0,.25,.5,.25)+.5,y=c(.75,.5,.75,1)))
w3 <- owin(poly=list(x=c(0,.25,.5,.25),y=c(.75,.5,.75,1)-.5))
w4 <- owin(poly=list(x=c(0,.25,.5,.25)+.5,y=c(.75,.5,.75,1)-.5))
w5 <- owin(poly=list(x=c(0.25,.5,.75,.5),y=c(.5,.25,.5,.75)))
wWrong <- owin(poly=list(x=c(0,1,1,0),y=c(.75,.5,.75,1)))
wList = lapply(seq_len(nfov), function(x){
    list("w1" = w1, "w2" = w2, "w3" = w3, "w4" = w4, "w5" = w5)
})
wList2 = lapply(seq_len(nfov), function(x){
    list("w1" = w1, "w2" = w2, "w3" = w3, "w4" = w4, "w5" = w5, "wWrong" = wWrong)
})
names(wList) = names(wList2) = seq_len(nfov)
hypFrame <- buildHyperFrame(df, coordVars = c("x", "y"), designVar = "fov")
hypFrame2 = addCell(hypFrame, wList)
#Register the parallel backend
nCores = 2
register(MulticoreParam(nCores))
test_check("spatrans")
