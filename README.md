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

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

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
    ## SmBIRDa    0.2636815 0.006749991 4.887131e-24 3.909705e-22
    ## SmAUX1a    0.3407263 0.005450551 5.694075e-22 2.277630e-20
    ## SmCYCA1;1a 0.3482799 0.006797150 6.187973e-19 1.650126e-17
    ## SmPFA2b    0.2278342 0.010168430 3.895420e-17 7.790840e-16
    ## SmCYCD3;3a 0.3903937 0.005992172 9.671855e-17 1.547497e-15
    ## SmSGNb     0.3397532 0.010935983 8.886705e-14 1.184894e-12

Finally write the results to a spreadsheet:

``` r
writeToXlsx(allModsNN, file = "myfile.xlsx")
```
