SMOPPIX: Single-MOlecule sPatial omics data analysed using the
Probabilistic IndeX
================
Stijn Hawinkel and Steven Maere

This repo provides code for analysing spatial transcriptomics and other
omics data on the single-molecule level using probabilistic indices. A
simple use-case is shown below, more extensive documentation can be
found in the vignette

The package can be installed from GitHub as follows:

``` r
library(devtools)
install_github("sthawinke/smoppix")
```

Load the package

``` r
library(smoppix)
```

Load a dataset, and convert it to a *spatstat* hyperframe.

``` r
data(Yang)
hypYang <- buildHyperFrame(Yang, coordVars = c("x", "y"), 
                           imageVars = c("day", "root", "section"))
```

    ## Found 29 unique images

Estimate the univariate nearest neighbour probabilistic index:

``` r
nnObj <- estPis(hypYang, pis = "nn", null = "background", verbose = FALSE)
```

Add a variance weighting function and plot it

``` r
nnObj <- addWeightFunction(nnObj, designVars = c("day", "root"))
plotWf(nnObj, pi = "nn")
```

![](README_files/figure-gfm/wf-1.png)<!-- -->

The inference step: fit linear mixed models, and show the most
significant results:

``` r
allModsNN <- fitLMMs(nnObj, fixedVars = "day", randomVars = "root")
```

    ## Fitted formula for pi nn:
    ## pi - 0.5 ~ 1 + day + (1 | root)

``` r
head(getResults(allModsNN, "nn", "Intercept"))
```

    ##             Estimate          SE         pVal         pAdj
    ## SmBIRDa    0.2634783 0.006729239 4.403772e-24 3.523017e-22
    ## SmAUX1a    0.3405152 0.005414905 4.630934e-22 1.852374e-20
    ## SmCYCA1;1a 0.3482347 0.006806002 6.349394e-19 1.693172e-17
    ## SmPFA2b    0.2273036 0.010080643 3.167296e-17 6.334592e-16
    ## SmCYCD3;3a 0.3903064 0.005979070 8.974696e-17 1.435951e-15
    ## SmSGNb     0.3399488 0.010944961 9.304194e-14 1.240559e-12

Finally write the results to a spreadsheet:

``` r
writeToXlsx(allModsNN, file = "myfile.xlsx")
```
