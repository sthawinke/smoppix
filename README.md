<img src='inst/Smoppix.jpg' align='centre' height='15%' width='15%'/>SMOPPIX:
Single-MOlecule sPatial omics data analysed using the Probabilistic
IndeX
================

This repo provides code for analyzing single-molecule spatial omics data
and cell type location data using probabilistic indices as introduced in
our [preprint](https://doi.org/10.1101/2025.05.20.654270). A simple
use-case is shown below, more extensive documentation can be found in
the vignette.

The package can be installed from Bioconductor as follows:

``` r
library(BiocManager)
install("smoppix")
```

The latest version can be installed from this Github repo as:

``` r
library(devtools)
install_github("sthawinke/smoppix", build_vignettes = TRUE)
```

Once installed, you can load the package

``` r
library(smoppix)
```

Patterns that can be detected by *smoppix*:

<img src='inst/CSR.jpg' align='centre' height='55%' width='55%'/>
<img src='inst/CSRbi.jpg' align='centre' height='55%' width='55%'/>

For illustration, we now load an example dataset, contained in the
package. It is in table format, so we first convert it to a *spatstat*
hyperframe.

``` r
data(Yang)
hypYang <- buildHyperFrame(Yang,
    coordVars = c("x", "y"),
    imageVars = c("day", "root", "section")
)
```

    ## Found 29 unique images

The number of unique images found is printed, make sure that this is
what you expected. Now to make completely sure the software has
understood us, we make an exploratory plot:

``` r
plotExplore(hypYang)
```

![](README_files/figure-gfm/explPlot-1.png)<!-- -->

As an example analysis, we estimate the univariate nearest neighbour
probabilistic index as a measure of aggregation (clustering):

``` r
nnObj <- estPis(hypYang, pis = "nn", null = "background", verbose = FALSE)
```

We add a variance weighting function to prepare for analysis and plot
it:

``` r
nnObj <- addWeightFunction(nnObj, designVars = c("day", "root"))
plotWf(nnObj, pi = "nn")
```

![](README_files/figure-gfm/wf-1.png)<!-- -->

As expected, the weight allotted to the point pattern increases with the
number of molecules in it: denser patterns provide more precise
estimates of localization patterns.

Next is the inference step: we fit linear mixed models, and show the
most significant results:

``` r
allModsNN <- fitLMMs(nnObj, fixedVars = "day", randomVars = "root")
```

    ## Fitted formula for pi nn:
    ## pi - 0.5 ~ 1 + day + (1 | root)

``` r
head(getResults(allModsNN, "nn", "Intercept"))
```

    ##             Estimate          SE         pVal         pAdj
    ## SmBIRDa    0.2631182 0.006767208 4.908700e-24 4.074221e-22
    ## SmAUX1a    0.3398707 0.005336197 2.837610e-22 1.177608e-20
    ## SmCYCA1;1a 0.3464801 0.006739431 3.671752e-19 1.015851e-17
    ## SmCYB2;4   0.3213716 0.008885568 8.940561e-18 1.855166e-16
    ## SmPFA2b    0.2270041 0.009784726 1.733721e-17 2.877977e-16
    ## SmCYCD3;3a 0.3905612 0.005936644 7.958008e-17 1.100858e-15

Let’s make it visual and plot the most significantly aggregated
transcripts:

``` r
plotTopResults(hypYang, allModsNN, pi = "nn")
```

![](README_files/figure-gfm/plotTopRes-1.png)<!-- -->

Finally write the results to a spreadsheet:

``` r
writeToXlsx(allModsNN, file = "myfile.xlsx")
```

A more extensive description of the *smoppix* functionality can be found
in the vignette, which can be accessed by calling

``` r
browseVignettes("smoppix")
```

A schematic representation of the *smoppix* pipeline is shown below:

<img src='inst/tikzPic.jpg' align='centre' height='100%' width='100%'/>
