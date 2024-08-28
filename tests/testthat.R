library(testthat)
library(smoppix)
library(spatstat.random)
library(BiocParallel)
n <- 5000 # number of molecules
ng <- 8 # number of genes
nfov <- 3 # Number of fields of view
conditions <- 4
# sample xy-coordinates in [0, 1]
x <- runif(n)
y <- runif(n)
# assign each molecule to some gene-cell pair
gs <- paste0("gene", seq(ng))
gene <- sample(gs, n, TRUE)
fov <- as.character(sample(nfov, n, TRUE))
condition <- as.character(sample(conditions, n, TRUE))
f = paste(fov, condition, sep = "_")
age = unsplit(lapply(split(integer(n), f = f), function(x) {rep(runif(1, 18, 98))}), f = f)
# construct data.frame of molecule coordinates
df <- data.frame(gene, x, y, fov, condition = condition, age = age)
# A list of point patterns
listPPP <- tapply(seq(nrow(df)), df$fov, function(i) {
    ppp(x = df$x[i], y = df$y[i], marks = df[i, c("gene", "condition", "fov", "age"), drop = FALSE])
}, simplify = FALSE)
# Regions of interest (roi): Diamond in the center plus four triangles
w1 <- owin(poly = list(x = c(0, 0.5, 1, 0.5), y = c(0.5, 0, 0.5, 1)))
w2 <- owin(poly = list(x = c(0, 0, 0.5), y = c(0.5, 0, 0)))
w3 <- owin(poly = list(x = c(0, 0, 0.5), y = c(1, 0.5, 1)))
w4 <- owin(poly = list(x = c(1, 1, 0.5), y = c(0.5, 1, 1)))
w5 <- owin(poly = list(x = c(1, 1, 0.5), y = c(0, 0.5, 0)))
wWrong <- owin(poly = list(x = c(0, 1, 1, 0), y = c(0.75, 0.5, 0.75, 1)))
hypFrame <- buildHyperFrame(df, coordVars = c("x", "y"), imageVars = c(
    "condition", "fov", "age"
), imageIdentifier = c("fov", "condition"))
nDesignFactors <- length(unique(hypFrame$image))
wList <- lapply(seq_len(nDesignFactors), function(x) {
    Li <- list(w1 = w1, w2 = w2, w3 = w3, w4 = w4, w5 = w5)
    names(Li) <- paste0(names(Li), "_", x)
    Li
})
wList2 <- lapply(seq_len(nDesignFactors), function(x) {
    list(w1 = w1, w2 = w2, w3 = w3, w4 = w4, w5 = w5, wWrong = wWrong)
})
unCells <- unlist(lapply(wList, names))
cellTypesDf <- data.frame(cell = unCells, cellType = sample(paste0("CellType_", LETTERS[seq_len(5)]),
    length(unCells),
    replace = TRUE
))
names(wList) <- names(wList2) <- rownames(hypFrame)
hypFrame2 <- addCell(hypFrame, wList, cellTypes = cellTypesDf, verbose = FALSE)
# The nuclei
n1 <- owin(poly = list(x = c(0.2, .4, 0.8, .4), y = c(.4, .2, .4, .8)))
n2 <- owin(poly = list(x = c(0.1, 0.1, .4), y = c(.4, .1, .1)))
n3 <- owin(poly = list(x = c(0.1, 0.1, .4), y = c(1, .75, 1)))
n4 <- owin(poly = list(x = c(1, 1, .6), y = c(.7, .9, .9)))
n5 <- owin(poly = list(x = c(.95, .95, .7), y = c(.1, .4, .1)))
nLarge <- owin(poly = list(x = c(1.1, 1.1, .7), y = c(.1, .4, .1)))
nList <- lapply(seq_len(nDesignFactors), function(x) {
    list("w1" = n1, "w2" = n2, "w3" = n3, "w4" = n4, "w5" = n5)
})
nList2 <- lapply(seq_len(nDesignFactors), function(x) {
  list("w1" = n1, "w2" = n2, "w3" = n3, "w4" = n4, "w5" = nLarge)
})
names(nList) <- names(nList2) <- rownames(hypFrame) # Matching names is necessary
hypFrame2 <- addNuclei(hypFrame2, nList, verbose = FALSE)
# Register the parallel backend
nCores <- 2
register(MulticoreParam(nCores))
# register(SerialParam()) # Switch on when mapping test coverage
pis <- c("nn", "nnPair", "edge", "centroid", "nnCell", "nnPairCell")
piEstsBG <- estPis(hypFrame2, pis = pis, null = "background", verbose = FALSE)
piEstsCSR <- estPis(hypFrame2, pis = pis, null = "CSR", verbose = FALSE)
piEstsBG2 <- estPis(hypFrame2[, c("ppp", "image", "tabObs")],
    pis = "nn", null = "background",
    verbose = FALSE
)
# Add weight functions
objBG <- addWeightFunction(piEstsBG, designVars = "condition")
objCSR <- addWeightFunction(piEstsCSR, designVars = "condition")
# Fit Yang models too
data(Yang)
hypYang <- buildHyperFrame(Yang, coordVars = c("x", "y"), imageVars = c(
    "day", "root", "section"
))
yangPims <- estPis(hypYang, features = getFeatures(hypYang)[12:18], 
                   pis = c("nn", "nnPair"), verbose = FALSE, nPointsAll = 4000)
yangPims <- addWeightFunction(yangPims, lowestLevelVar = "section")
data(Eng)
hypEng <- buildHyperFrame(Eng, coordVars = c("x", "y"), imageVars = c("fov", "experiment"))
hypEng <- addCell(hypEng, EngRois, verbose = FALSE)
test_check("smoppix")
