[![DOI](https://zenodo.org/badge/75324979.svg)](https://zenodo.org/badge/latestdoi/75324979)

# Effective Degrees of Freedom of the Pearson's Correlation Coefficient under Serial Correlation

## Highlights
* Autocorrelation biases the standard error of Pearson's correlation and breaks the variance-stabilising property of Fisher's transformation.
* Severity of resting state fMRI autocorrelation varies systematically with region of interest size, and is heterogeneous over subjects.
* Commonly used methods (see `mis` directory) to adjust correlation standard errors are themselves biased when true correlation is non-zero due to a confounding effect.
* We propose a “xDF” method to provide accurate estimates of the variance of Pearson’s correlation -- before or after Fisher’s transformation -- that considers auto-correlation of each time series as well as instantaneous and lagged cross-correlation.
* Accounting for the autocorrelation in resting-state functional connectivity considerably alters the graph theoretical description of human connectome.



## Table of contents
* [Introduction](#introduction)
* [Configurations](#Configurations)
  * [Dependencies](##Dependencies)
  * [Octave](##Octave)
* [Examples](#Examples)
  * [Using xDF](#xxDF)
  * [Constructing Functional Connectivity (FC) Maps](#FC)
    * [FDR-based Statistically Thresholded Functional Connectivity](#STFC)
    * [CE-based proportionally Thresholded Functional Connectivity](#CEFC)
    * [Unthresholded Functional Connectivity](#UFC)
* [Simulating time series of arbitrary correlation and autocorrelation structure](#Sim)
  * [Correlated but White Time Series](#CW)
  * [Uncorrelated but Autocorrelated Time Series](#UA)
  * [Correlated and Autocorrelated Time Series](#CA)




## Introduction <a name="introduction"></a>
Collection of scripts to implement the xDF method introduced in

*Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. "Effective Degrees of Freedom of the Pearson's Correlation Coefficient under Serial Correlation." bioRxiv (2018): 453795.*

The `xDF.m` may be used to:
* Estimate the variance of Pearson's correlation coefficients
* Calculate the z-statistics maps of functional connectivity
* Estimate the p-values for such correlation coefficients


## Configurations <a name="configurations"></a>
For now, the xDF has only been implemented in MATLAB. Although we will be releasing the Python version in a near future. You need MATLAB statistical toolbox to run the script.

To clone the repository use:
```
git clone https://github.com/asoroosh/xDF.git
```
or alternatively download the [zip file](https://github.com/asoroosh/xDF/archive/master.zip).

### Dependencies <a name="dependencies"></a>

The `xDF.m` should work without requiring any external function. However it is comprised of several internal modules which are also available individually in `mis/`:

* `xC_fft.m`: approximates the cross correlations using Wiener–Khinchin theorem
* `AC_fft.m`: approximates the autocorrelations using Wiener–Khinchin theorem


### Octave <a name="octave"></a>
xDF is also available for Octave via `xDF_octave.m`. Note that you require statistics package to run the script in Octave:
```
pkg install -forge io
pkg install -forge statistics
```
We have only tested the script on Octave 4.4.1.


## Examples <a name="examples"></a>

### Using xDF <a name="xxDF"></a>
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

### Constructing Functional Connectivity (FC) Maps <a name="fc"></a>
Here we show how you can use xDF to form functional connectivity of BOLD signals of `I` region of interest and `T` time-points. The examples cover both statistically and proportionally thresholded FCs as well as unthresholded FCs.
#### FDR-based Statistically Thresholded Functional Connectivity <a name="stfc"></a>

```
[VarRho,Stat]=xDF(ts,T,'truncate','adaptive','TVOff')
FC = fdr_bh(Stat.p).*Stat.z; %FC of IxI
```
Function `fdr_bh` is an external function [[+]](https://uk.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh?focused=5807896&tab=function). It can also be found in `.../xDF/FCThresholding/StatThresholding/`
#### CE-based proportionally Thresholded Functional Connectivity <a name="cefc"></a>
```
[VarRho,Stat]=xDF(ts,T,'truncate','adaptive','TVOff')
densrng = 0.01:0.01:0.50;
[~,CE_den]=CostEff_bin(Stat.z,densrng)
FC = threshold_proportional(Stat.z,CE_den); %FC of IxI
```

Function `threshold_proportional` is an external function from Brain Connectivity Toolbox [[+]](https://sites.google.com/site/bctnet/Home/functions).

#### Unthresholded Functional Connectivity <a name="ufc"></a>
```
[VarRho,Stat]=xDF(ts,T,'truncate','adaptive','TVOff')
FC = Stat.z; %FC of IxI
```

## Simulating time series of arbitrary correlation and autocorrelation structure <a name="sim"></a>

If you are interested in reproducing results in the paper or sanity check the xDF. You can simulate N time series of desired correlation matrix of `C` and autocorrelation of `A` using function `corrautocorr`. We are showing example of three scenarios for pair time series of length `T`.

### Correlated but White Time Series <a name="cw"></a>

```
>> ts = corrautocorr([0 0],0.9,eye(T),T);
>> corr(ts')

ans =

    1.0000    0.9014
    0.9014    1.0000
```

### Uncorrelated but Autocorrelated Time Series <a name="ua"></a>

```
C_T1 = MakeMeCovMat([0.9:-.1:.1],T); %autocorrelation matrix of first time series (AR9=0.9:0.1)
C_T2 = MakeMeCovMat(0.4,T); %autocorrelation matrix of first time series (AR1 = 0.4)
A = cat(3,C_T1,C_T2); % TxTxI time series autocorrelation matrices
>> ts = corrautocorr([0 0],0,A,T);
>> rhohat = corr(ts')

rhohat =

    1.0000   -0.0350
   -0.0350    1.0000

AC_ts = AC_fft(ts,T); % estimate the autocorrelations

ac_x  = AC_ts(1,1:T-1);
ac_y  = AC_ts(2,1:T-1);

ac_x =

    1.0000    0.8962    0.7922    0.6924    0.6069    0.5357    0.4651    0.3770    0.2810    

ac_y =

    1.0000    0.4138    0.0143    0.0269    0.0765    0.0509    0.0116   -0.0213   -0.0587   
```

### Correlated and Autocorrelated Time Series <a name="ca"></a>

```
C_T1 = MakeMeCovMat([0.9:-.1:.1],T); %autocorrelation matrix of first time series (AR9=0.9:0.1)
C_T2 = MakeMeCovMat(0.4,T); %autocorrelation matrix of second time series (AR1 = 0.4)
rho = 0.4;
A = cat(3,C_T1,C_T2); % TxTxI time series autocorrelation matrices
C = [1 rho; rho 1]; % IxI correlation matrix

ts = corrautocorr([0 0],C,A,T); %simulates the time series using Cholesky decomposition

AC_ts = AC_fft(ts,T); % estimate the autocorrelations

ac_x  = AC_ts(1,1:T-1);
ac_y  = AC_ts(2,1:T-1);

ac_x =

  1.0000    0.9021    0.8049    0.6960    0.5850    0.4955    0.4091    0.3137    0.2109

ac_y =

  1.0000    0.4885    0.0580   -0.0824   -0.1152   -0.1112   -0.0734   -0.0893   -0.0139

rhohat = corr(ts(1,:)',ts(2,:)')

rhohat =

  0.1239

```

It is important to note that when two time series are highly correlated _and_ they have different autocorrelation structure, the sample correlation is underestimated. Importantly, this shouldn't be confused with the fact that the sample correlations are approximately unbiased even if the time series are autocorrelated. This could be corrected using Eq. S9 in Section S3.1 as following

```
Sigma_X  = toeplitz(ac_x);
Sigma_Y  = toeplitz(ac_y);

Kx = chol(Sigma_X);
Ky = chol(Sigma_Y);
rhohat_corrected = rhohat./(trace(Kx*Ky')./T);

rhohat_corrected =

  0.3888
```
Obviously, these are just quick examples; more accurate estimates should be achieved via hundreds of iterations.
