
# 0.1.0

- Published on github

# 0.1.1

- Use regular rather than squared distances to calculate nearest
  neighbour PIs

# 0.99.0

- Submission to BioConductor

# 0.99.4

- Prediction of cell types

# 0.99.9

- Prediction of cell types

# 0.99.9

- Bugfixes in hyperframe building

# 0.99.10

- Reducing computation time of the examples

# 0.99.17

- Bugfix in calcnnPi

# 0.99.19

- Shorter example times

# 0.99.20

- Reducing running times of examples to pass R CMD CHECK

# 0.99.21

- Bugfix in buildDataFrame

# 0.99.23

- Added addNuclei() function for adding nuclei as spatstat::owin

# 0.99.24

- Nuclei can be plotted using plotExplore()

# 0.99.25

- ggplot2 objects are now returned for user customization
- Repeated return() calls are omitted
- NEWS file was updated to reflect all changes

# 0.99.26

- Version bump with serial processing to find bug occurring only on
  server

# 0.99.30

- Reduce computation time of tests and vignettes in view of BioConductor
  server warnings
- Automate choice of parallel backend for different operating systems in
  vignette

# 0.99.36

- Simplify documentation and reducing example runtimes by pooling
  manpages using the rdname commend

# 0.99.39

- Further round of reducing documentation and test time by pooling
  examples
- Test windows on a single core

# 0.99.40

- Use DOIs in URLs following winbuilder
- Use makeCluster() on windows in the examples too

# 0.99.41

- Clean up manpages for functions sharing the same example set

# 0.99.42

- Reviewed all manpages, clarifying user options

# 0.99.43

- Fixed bug on windows using bpnworkers()

# 0.99.44

- Limit number of features to avoid windows timeout

# 1.1.1

- Allow to pass own colours to plotExplore

# 1.1.2

- Calculate nnPair PIs only between desired features

# 1.1.3

- Some speedups in the buildDataFrame function
- With reference to bioRxiv preprint

# 1.1.4

- Optimize addWeightFunction

# 1.1.5

- Export loadBalanceBplappy for use in other packages
- Mention Bioconductor installer in README
