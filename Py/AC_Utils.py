#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 16:53:44 2018

@author: sorooshafyouni
University of Oxford, 2019
srafyouni@gmail.com
"""

def AC_fft(Y,T):
    import numpy as np 
    
    
    if np.shape(Y)[1]!=T:
        print('AC_fft::: Input should be in IxT form, the matrix was transposed.')
        Y = np.transpose(Y)
        
    print("AC_fft::: Demean along T")
    mY2 = np.mean(Y,axis=1)
    Y = Y-np.transpose(np.tile(mY2,(T,1)));
    
    nfft    = nextpow2(2*T-1); #zero-pad the hell out!
    yfft    = np.fft.fft(Y,n=nfft,axis=1); #be careful with the dimensions
    ACOV = np.real(np.fft.ifft(yfft*np.conj(yfft),axis=1));
    ACOV = ACOV[:,0:T-1];
    
    Norm = np.sum(np.abs(Y)**2,axis=1);
    Norm = np.transpose(np.tile(Norm,(T-1,1)));
    xAC  = ACOV/Norm; #normalise the COVs   
    
    bnd  = (np.sqrt(2)*1.3859)/np.sqrt(T); #assumes normality for AC
    CI   = [-bnd,bnd];
    
    return xAC,CI


def xC_fft(Y,T, mxL = [], **kwargs):
    import numpy as np

# ***********************************
# This should be checked! There shouldn't be any complex numbers!!
#__main__:74: ComplexWarning: Casting complex values to real discards the imaginary part
# ***********************************
    
    if np.shape(Y)[1]!=T:
        print('xC_fft::: Input should be in IxT form, the matrix was transposed.')
        Y = np.transpose(Y)
    
    if not np.size(mxL): 
        mxL = T;
    
    I = np.shape(Y)[0]
    
    print("AC_fft::: Demean along T")
    mY2 = np.mean(Y,axis=1)
    Y   = Y-np.transpose(np.tile(mY2,(T,1)));    
    
    nfft    = nextpow2(2*T-1); #zero-pad the hell out!
    yfft    = np.fft.fft(Y,n=nfft,axis=1); #be careful with the dimensions    

    mxLcc = (mxL-1)*2+1
    xC    = np.empty([I,I,mxLcc]);
        
    XX = np.where(np.triu(np.ones([I,I]),1))[0]
    YY = np.where(np.triu(np.ones([I,I]),1))[1]
        
    for i in np.arange(0,np.size(XX)): #loop around edges
        
        xC0 = np.fft.ifft(yfft[XX[i],:]*np.conj(yfft[YY[i],:]))
        xC0 = np.concatenate((xC0[-mxL+1:],xC0[0:mxL]))                
        xC0 = np.fliplr([xC0])[0]        
        Norm = np.sqrt(np.sum(np.abs(Y[XX[i],:])**2)*np.sum(np.abs(Y[YY[i],:])**2))     
                
        xC0 = xC0/Norm;        
        xC[XX[i],YY[i],:]   = xC0
        del xC0
    
    xC = xC + np.transpose(xC,(1,0,2));
    
    lidx = np.arange(-(mxL-1),mxL);
    
    return xC,lidx

def nextpow2(x):
    """
    nextpow2 Next higher power of 2.
    nextpow2(N) returns the first P such that P >= abs(N).  It is
    often useful for finding the nearest power of two sequence
    length for FFT operations.
    """
    import math
    return 1 if x == 0 else 2**math.ceil(math.log2(x))

#import scipy.io
#import numpy as np
#V = '/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_124422_OnlyMTS.mat'
#mat = scipy.io.loadmat(V)
#mts = mat['mts']
#T = np.shape(mts)[0]
#[AC,CI] = AC_fft(mts,T)
#[xC,lid] = xC_fft(mts,T)

