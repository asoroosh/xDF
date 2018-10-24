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

### Using xDF
Suppose `X` is a matrix of size `IxT` where `I` is number of regions and `T` is number of data-points,

1) Estimating the variance *without* any regularisation on the autocorrelation function:

```
[V,Stat]=xDF(X,T)
VarRho % is matrix of IxI where each element is variance of correlation coefficient between corresponding elements.
Stat.p % matrix of IxI of uncorrected p-values
Stat.z % matrix of IxI of z-scores
```
2) Estimating the variance with Tukey tapering method, on cut off `M=sqrt(T)` as suggested in Chatfield 2016.

`[VarRho,Stat]=xDF(X,T,'taper','tukey',sqrt(T))`

Similar regularisation with `M=2*sqrt(T)` as in Woolrich et al 2001:

`[VarRho,Stat]=xDF(X,T,'taper','tukey',2*sqrt(T))`

3) Estimating the variance with adaptive truncation method. Which we showed that generates the most accurate estimates.

`[VarRho,Stat]=xDF(X,T,'truncate','adaptive','verbose')`

4) Estimating the variance with shrinking as tapering method & without controlling the variance.

`[VarRho,Stat]=xDF(X,T,'truncate','adaptive','TVOff')`

For more options see the usage section of the function.

### Constructing Functional Connectivity (FC) Maps
Here we show how you can use xDF to form functional connectivity of BOLD signals of `I` region of interest and `T` time-points. The examples cover both statistically and proportionally thresholded FCs as well as unthresholded FCs.
#### FDR-based Statistically Thresholded Functional Connectivity

```
[VarRho,Stat]=xDF(ts,T,'truncate','adaptive','TVOff')
FC = fdr_bh(Stat.p).*Stat.z; %FC of IxI
```
Function `fdr_bh` is an external function [[+]](https://uk.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh?focused=5807896&tab=function). It can also be found in `.../xDF/FCThresholding/StatThresholding/`
#### CE-based proportionally Thresholded Functional Connectivity
```
[VarRho,Stat]=xDF(ts,T,'truncate','adaptive','TVOff')
densrng = 0.01:0.01:0.50;
[~,CE_den]=CostEff_bin(Stat.z,densrng)
FC = threshold_proportional(Stat.z,CE_den); %FC of IxI
```

Function `threshold_proportional` is an external function from Brain Connectivity Toolbox [[+]](https://sites.google.com/site/bctnet/Home/functions).

#### Unthresholded Functional Connectivity
```
[VarRho,Stat]=xDF(ts,T,'truncate','adaptive','TVOff')
FC = Stat.z; %FC of IxI
```
