# xDF

# Table of contents
1. [Introduction](#introduction)
2. [Configurations](#Configurations)
  1. [Dependencies](#Dependencies)
3. [Examples](#Examples)
  1. [Using xDF](#xDF)
  2. [Constructing Functional Connectivity (FC) Maps](#FC)
    1. [FDR-based Statistically Thresholded Functional Connectivity](#STFC)
    2. [CE-based proportionally Thresholded Functional Connectivity](#CEFC)
    3. [Unthresholded Functional Connectivity](#UFC)
4. [Simulating time series of arbitrary correlation and autocorrelation structure](#Sim)
  1. [Correlated but White Time Series](#CW)
  2. [Uncorrelated but Autocorrelated Time Series](#UA)
  3. [Correlated and Autocorrelated Time Series](#CA)




## Introduction <a name="introduction"></a>
Collection of scripts to implement the xDF method introduced in

*Afyouni, Smith & Nichols, Effective Degrees of Freedom of the Pearson's Correlation Coefficient under Serial Correlation, BioRxiv DOI: XX - XX - XX*

The `xDF.m` may be used to:
* Estimate the variance of Pearson's correlation coefficients
* Calculate the z-statistics maps of functional connectivity
* Estimate the p-values for such correlation coefficients


## Configurations <a name="configurations"></a>
For now, the xDF has only been implemented in MATLAB. Although we will be releasing the Python version in a near future. You need MATLAB statistical toolbox to run the script.

### Dependencies <a name="dependencies"></a>

The `xDF.m` should work without requiring any external function. However it was comprised of several modules which are available individually in `Aux/`

* `xC_fft.m`: approximates the cross correlations using Wiener–Khinchin theorem
* `AC_fft.m`: approximates the autocorrelations using Wiener–Khinchin theorem

## Examples <a name="Examples"></a>

### Using xDF <a name="xDF"></a>
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

### Constructing Functional Connectivity (FC) Maps <a name="FC"></a>
Here we show how you can use xDF to form functional connectivity of BOLD signals of `I` region of interest and `T` time-points. The examples cover both statistically and proportionally thresholded FCs as well as unthresholded FCs.
#### FDR-based Statistically Thresholded Functional Connectivity <a name="STFC"></a>

```
[VarRho,Stat]=xDF(ts,T,'truncate','adaptive','TVOff')
FC = fdr_bh(Stat.p).*Stat.z; %FC of IxI
```
Function `fdr_bh` is an external function [[+]](https://uk.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh?focused=5807896&tab=function). It can also be found in `.../xDF/FCThresholding/StatThresholding/`
#### CE-based proportionally Thresholded Functional Connectivity <a name="CEFC"></a>
```
[VarRho,Stat]=xDF(ts,T,'truncate','adaptive','TVOff')
densrng = 0.01:0.01:0.50;
[~,CE_den]=CostEff_bin(Stat.z,densrng)
FC = threshold_proportional(Stat.z,CE_den); %FC of IxI
```

Function `threshold_proportional` is an external function from Brain Connectivity Toolbox [[+]](https://sites.google.com/site/bctnet/Home/functions).

#### Unthresholded Functional Connectivity <a name="UFC"></a>
```
[VarRho,Stat]=xDF(ts,T,'truncate','adaptive','TVOff')
FC = Stat.z; %FC of IxI
```

## Simulating time series of arbitrary correlation and autocorrelation structure <a name="Sim"></a>

If you are interested in reproducing results in the paper or sanity check the xDF. You can simulate N time series of desired correlation matrix of `C` and autocorrelation of `A` using function `corrautocorr`. We are showing example of three scenarios for pair time series of length `T`.

### Correlated but White Time Series <a name="CW"></a>

```
>> ts = corrautocorr([0 0],0.9,eye(T),T);
>> corr(ts')

ans =

    1.0000    0.9014
    0.9014    1.0000
```

### Uncorrelated but Autocorrelated Time Series <a name="UA"></a>

```
C_T1 = MakeMeCovMat([0.9:-.1:.1],T); %autocorrelation matrix of first time series (AR9=0.9:0.1)
C_T2 = MakeMeCovMat(0.4,T); %autocorrelation matrix of first time series (AR1 = 0.4)
A = cat(3,C_T1,C_T2); % TxTxI time series autocorrelation matrices
>> ts = corrautocorr([0 0],0,A,T);
>> corr(ts')

ans =

    1.0000   -0.0350
   -0.0350    1.0000

>> AC1 = autocorr(ts(1,:),9)

AC1 =

    1.0000    0.8962    0.7922    0.6924    0.6069    0.5357    0.4651    0.3770    0.2810    0.1834

>> AC2 = autocorr(ts(2,:),9)

AC2 =

    1.0000    0.4138    0.0143    0.0269    0.0765    0.0509    0.0116   -0.0213   -0.0587   -0.0593
```

### Correlated and Autocorrelated Time Series <a name="CA"></a>

```
C_T1 = MakeMeCovMat([0.9:-.1:.1],T); %autocorrelation matrix of first time series (AR9=0.9:0.1)
C_T2 = MakeMeCovMat(0.4,T); %autocorrelation matrix of first time series (AR1 = 0.4)
rho = 0.4;
A = cat(3,C_T1,C_T2); % TxTxI time series autocorrelation matrices
C = [1 rho; rho 1]; % IxI correlation matrix

ts = corrautocorr([0 0],C,A,T); %simulates the time series using Cholesky decomposition

>> corr(ts')

ans =

    1.0000    0.1938
    0.1938    1.0000

>> AC1 = autocorr(ts(1,:),9)

AC1 =

    1.0000    0.8992    0.8027    0.6958    0.5933    0.4895    0.3916    0.2960    0.1908    0.0816

>> AC2 = autocorr(ts(2,:),9)

AC2 =

    1.0000    0.4095   -0.0211   -0.0300    0.0033    0.0722    0.0520    0.0565    0.0432    0.0025
```
It is important to note that in presence of autocorrelation, the estimated correlation is usually lower than the true correlation. See Section S3.1 in the paper.
