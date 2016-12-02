fMRI BOLD signals are notorious in terms of violating the dependency presumption in statistical tests. This dependency result in larger variance and therefore may bias the statsitical hypothesis tests as well as analysis of topological features of the human connectome. This toolbox was designed to find a correction factor to scale the biased variance of correlation coefficients between pair fMRI BOLD signals. 

# HetBiv
HetBiv correction factors are a family of correctio methods which estimate the effective degree-of-freedom to address the inflation in variance due to dependency. The three main methods of this family differs in terms of dependency assumptions *between* two time series under study. 

## Hetrogenous Autocorrelation Effect
The effect of auctorrelation is not spatially homogenous across the grey matter. This hetrogeneity bolds the importance of using *local* correction factors. In other words, the correction factors are estimated for edge between each pair nodes in the functional connectivity.

## HetBiv for Independent (un-correlated) time series
In this case, the assumption is that the two time series are independent. To estimate the correction factor, in this case, the function *HetBiv_Short.m* or *HetBiv_Short_Fast.m* should be used. This function produces a correction factor (CF) which can later be used to estimate effective degree-of-freedom, EDoF=DoF/CF;

## HetBiv for 0-lag dependent (correlated) time series


## HetBiv for n-lag dependent (correlated & cross-correlated) time series
