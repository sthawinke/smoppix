library(testthat)
library(spatrans)
library(spatstat.random)
n <- 1e3 # number of molecules
ng <- 50 # number of genes
nc <- 20 # number of cells
nfov = 3 # Number of fields of view
# sample xy-coordinates in [0, 1]
x <- runif(n)
y <- runif(n)
# assign each molecule to some gene-cell pair
gs <- paste0("gene", seq(ng))
cs <- paste0("cell", seq(nc))
gene <- sample(gs, n, TRUE)
cell <- sample(cs, n, TRUE)
fov <- sample(nfov, n, TRUE)
# assure gene & cell are factors so that
# missing observations aren't dropped
gene <- factor(gene, gs)
cell <- factor(cell, cs)
# construct data.frame of molecule coordinates
df <- data.frame(gene, cell, x, y, fov)
#A list of point patterns
listPPP = tapply(seq(nrow(df)), df$fov, function(i){
    ppp(x = df$x[i], y = df$y[i], marks = df[i, c("gene", "cell")])
}, simplify = FALSE)
test_check("spatrans")
