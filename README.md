fMRI BOLD signals are notorious in terms of violating the dependency presumption in statistical tests. This dependency result in larger variance and therefore may bias the statsitical hypothesis tests as well as analysis of topological features of the human connectome. This toolbox was designed to find a correction factor to scale the biased variance of correlation coefficients between pair fMRI BOLD signals. 

# HetBiv
HetBiv correction factors are a family of correctio methods which estimate the effective degree-of-freedom to address the inflation in variance due to dependency. The three main methods of this family differs in terms of dependency assumptions *between* two time series under study. 

## Hetrogenous Autocorrelation Effect
The effect of auctorrelation is not spatially homogenous across the grey matter. This hetrogeneity bolds the importance of using *local* correction factors. In other words, the correction factors are estimated for edge between each pair nodes in the functional connectivity.

## HetBiv for Independent (un-correlated) time series
In this case, the assumption is that the two time series are independent. To estimate the correction factor, in this case, the function `HetBiv_Short.m` or `HetBiv_Short_Fast.m` should be used. This function produces a correction factor (CF) which can later be used to estimate effective degree-of-freedom, EDoF=DoF/CF;

## HetBiv for 0-lag dependent (correlated) time series
In this case, the assumption is that the 0 lag versions of the two time series are dependent and no cross-correlation exist between the two. To calculate the correction factor for this case, `HetBiv_ShortBC_Xior.m` and `HetBiv_ShortBC_Xior_Fast.m` should be used. Similarly, the EDoF is calculated as DoF/CF. 

## HetBiv for n-lag dependent (cross-correlated) time series
The last version and most sophisticated version of HetBiv correction factor assumes that the time series are not only dependent in their 0-lag versions, but also there are dependency between lagged versions of the two time series. However, in this case it is essential to use the *curbing factors* to avoid over-estimation of dependencies in far lags. The provisional validations suggest that the curbing factor depends on the length of the time series. However, they consistently suggest the smaller curbing factor performs better than the larger curbing factors. The correction factors can be estimated by `HetBiv_ShortBC.m` and `HetBiv_ShortBC_Fast.m`. 

## Usage
`CF=HetBiv_*(ts1,ts2,ndpr,howfar,figflag)`
`ts1 & ts2:` Two time series from node 1 and node 2. Note that the time series should have been already zero-meaned. 
`ndpr:` Number of time points. Note that the the number of time points should be equal in ts1 and ts2. 
`howfar:` Curbing factor which should be between 0 and 1. howfar=0 means 0lags and 1 indicated the whole time series. 
`figfalg:` Should be set to 1 if you are interested in toeplitz and covariance matrices, otherwise should be set to 0. 
