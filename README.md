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
    ## SmBIRDa    0.2637957 0.006708042 4.197715e-24 3.358172e-22
    ## SmAUX1a    0.3404423 0.005404623 4.353388e-22 1.741355e-20
    ## SmCYCA1;1a 0.3481950 0.006805181 6.287441e-19 1.676651e-17
    ## SmPFA2b    0.2274859 0.010084754 3.234472e-17 6.468944e-16
    ## SmCYCD3;3a 0.3901479 0.005949646 7.648429e-17 1.223749e-15
    ## SmSGNb     0.3401097 0.010919376 9.029209e-14 1.203895e-12

Finally write the results to a spreadsheet:

``` r
writeToXlsx(allModsNN, file = "myfile.xlsx")
```
