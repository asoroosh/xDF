# xDF

## Introduction
Collection of scripts to implement the xDF method introduced in

*Afyouni, Smith & Nichols, Effective Degrees of Freedom of the Pearson's Correlation Coefficient under Serial Correlation, BioRxiv DOI: XX - XX - XX*

The `xDF.m` may be used to:
* Estimate the variance of Pearson's correlation coefficients
* Calculate the z-statistics maps of functional connectivity
* Estimate the p-values for such correlation coefficients


## Configurations
For now, the xDF has only been implemented in MATLAB. Although we will be releasing the Python version in a near future. You need MATLAB statistical toolbox to run the script.

### Dependencies

The `xDF.m` should work without requiring any external function. However it was comprised of several modules which are available individually in `Aux/`

* `xC_fft.m`: approximates the cross correlations using Wiener–Khinchin theorem
* `AC_fft.m`: approximates the autocorrelations using Wiener–Khinchin theorem

## Examples

Suppose `X` is a matrix of size `IxT` where `I` is number of regions and `T` is number of data-points,

1) Estimating the variance *without* any tapering method

`[V,Stat]=xDF(X,T)`

2) Estimating the variance with Tukey tapering method

`[V,Stat]=xDF(X,T,'taper','tukey',sqrt(T))`

3) Estimating the variance with shrinking tapering method. We showed that `..,'taper','shrink',..` generates the most accurate estimates.

```[V,Stat]=xDF(X,T,'taper','shrink','verbose')```

4) Estimating the variance with shrinking as tapering method & without controlling the variance.

`[V,Stat]=xDF(X,T,'taper','shrink','TVOff')`

For more options see the usage section of the function. 
