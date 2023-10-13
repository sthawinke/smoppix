context("Test input spatrans package")
n <- 1e3  # number of molecules
ng <- 50  # number of genes
nc <- 20  # number of cells
nfov = 3#Number of fields of view
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
testthat("Reading in data proceeds without errors", {
  expect_message(hypFrame <- buildHyperFrame(df, coordVars = c("x", "y"), designVar = "fov"))
})
