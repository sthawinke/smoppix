SMOPPIX: Single-MOlecule sPatial omics data analysed using the
Probabilistic IndeX
================
Stijn Hawinkel

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

Load a dataset, convert it to a *spatstat* hyperframe, adn make an
explorative plot:

``` r
data(Yang)
hypYang <- buildHyperFrame(Yang, coordVars = c("x", "y"), 
                           imageVars = c("day", "root", "section"))
```

    ## Found 29 unique images

``` r
plotExplore(hypYang)
```

![](README_files/figure-gfm/loadYang-1.png)<!-- -->

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
    ## SmBIRDa    0.2638346 0.006757913 5.128128e-24 4.102502e-22
    ## SmAUX1a    0.3404722 0.005479555 6.275770e-22 2.510308e-20
    ## SmCYCA1;1a 0.3482532 0.006784812 5.879442e-19 1.567851e-17
    ## SmPFA2b    0.2276867 0.010201285 4.104503e-17 8.209007e-16
    ## SmCYCD3;3a 0.3904622 0.006009758 1.057280e-16 1.691648e-15
    ## SmSGNb     0.3398192 0.010948443 9.202101e-14 1.226947e-12

Finally write the results to a spreadsheet:

``` r
writeToXlsx(allModsNN, file = "myfile.xlsx")
```
