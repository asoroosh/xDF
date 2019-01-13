#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 13:31:32 2019

@author: sorooshafyouni
University of Oxford, 2019
"""
import numpy as np
from AC_Utils import *
from MatMan import *
import os,sys
import scipy.stats as sp

def xDF_Calc(ts,T,copy=True,verbose=False):
    
    if verbose: blockPrint()
    
    if copy: #Make sure you are not messing around with the original time series
        ts = ts.copy()
    
    if np.shape(ts)[1]!=T:
        print('xDF::: Input should be in IxT form, the matrix was transposed.')
        ts = np.transpose(ts)
    
    N  = np.shape(ts)[0];
    
    ts_std = np.std(ts,axis=1)
    ts     = ts/np.transpose(np.tile(ts_std,(T,1)));   #standardise

    #Corr----------------------------------------------------------------------
    rho    = CorrMat(ts,T)
    #Autocorr------------------------------------------------------------------
    [ac,CI] = AC_fft(ts,T); 
    ac   = ac[:,1:T-1]; #The last element of ACF is rubbish, the first one is 1, so why bother?!
    nLg  = T-2;  
    
    #Cross-corr---------------------------------------------------------------- 
    [xcf,lid] = xC_fft(ts,T);
    
    xc_p      = xcf[:,:,1:T-1];
    xc_p      = np.flip(xc_p,axis = 2); #positive-lag xcorrs
    xc_n      = xcf[:,:,T:-1];        #negative-lag xcorrs
    
    #Monster Equation---------------------------------------------------------
    wgt     = np.arange(nLg,0,-1)
    wgtm2   = np.tile((np.tile(wgt,[N,1])),[N,1]);
    wgtm3   = np.reshape(wgtm2,[N,N,np.size(wgt)]); #this is shit, eats all the memory!
    Tp      = T-1;
    
    print(np.shape(wgt))
    print(np.shape(wgtm2))
    print(np.shape(wgtm3))
    print(np.shape(ac))
    print(np.shape(xc_p))
    print(np.shape(xc_n))
    """
     VarHatRho = (Tp*(1-rho.^2).^2 ...
     +   rho.^2 .* sum(wgtm3 .* (SumMat(ac.^2,nLg)  +  xc_p.^2 + xc_n.^2),3)...         %1 2 4
     -   2.*rho .* sum(wgtm3 .* (SumMat(ac,nLg)    .* (xc_p    + xc_n))  ,3)...         % 5 6 7 8
     +   2      .* sum(wgtm3 .* (ProdMat(ac,nLg)    + (xc_p   .* xc_n))  ,3))./(T^2);   % 3 9 
    """
    
    # Da Equation!--------------------
    VarHatRho = (Tp*(1-rho**2)**2 \
                  +   rho**2  * np.sum(wgtm3 * (SumMat(ac**2,nLg)  +  xc_p**2 + xc_n**2),axis=2)\
                  -   2*rho   * np.sum(wgtm3 * (SumMat(ac,nLg)     * (xc_p    + xc_n))  ,axis=2)\
                  +   2       * np.sum(wgtm3 * (ProdMat(ac,nLg)    + (xc_p    * xc_n))  ,axis=2))/(T**2)    
 
    # diagonal is rubbish;
    VarHatRho[range(N),range(N)] = 0;
    
    #------- Test Stat-----------------------
    # Well, these are all Matlab and pretty useless -- copy pasted them just in case...
    #Pearson's turf -- We don't really wanna go there, eh?
    #rz      = rho./sqrt((ASAt));     %abs(ASAt), because it is possible to get negative ASAt!
    #r_pval  = 2 * normcdf(-abs(rz)); %both tails
    #r_pval(1:nn+1:end) = 0;          %NaN screws up everything, so get rid of the diag, but becareful here. 
    
    #Our turf--------------------------------
    rf      = np.arctanh(rho)
    sf      = VarHatRho/((1-rho**2)**2)    #delta method; make sure the N is correct! So they cancel out.
    rzf     = rf/np.sqrt(sf)
    rzf[range(N),range(N)] = 0;
    f_pval  = 2 * sp.norm.cdf(-abs(rzf));  #both tails
    f_pval[range(N),range(N)] = 0;         #NaN screws up everything, so get rid of the diag, but becareful here. 

    





    if verbose: enablePrint()
##############################################################################
##############################################################################
##############################################################################
    
# Disable verbose
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore verbose
def enablePrint():
    sys.stdout = sys.__stdout__    
################ TESTING PART:    
import scipy.io
#import numpy as np
V = '/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_124422_OnlyMTS.mat'
mat = scipy.io.loadmat(V)   
mts = mat['mts']
T = 1200
xDF_Calc(mts,T)