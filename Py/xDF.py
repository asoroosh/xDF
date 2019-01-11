#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 13:31:32 2019

@author: sorooshafyouni
University of Oxford, 2019
"""

def xDF_Calc(ts,T):
    import numpy as np
    import AC_Utils as acu
    
    if np.shape(ts)[1]!=T:
        print('xDF::: Input should be in IxT form, the matrix was transposed.')
        ts = np.transpose(ts)
    
    
    N  = np.shape(ts)[0];
    
    ts_std = np.std(ts,axis=1)
    ts     = ts/np.transpose(np.tile(ts_std,(T,1)));   #standardise

    #Corr----------------------------------------------------------------------
    rho    = CorrMats(ts)
    #Autocorr------------------------------------------------------------------
    [ac,CI] = acu.AC_fft(ts,T); 
    ac   = ac[:,1:T-2]; #The last element of ACF is rubbish, the first one is 1, so why bother?!
    nLg  = T-2;  
       
    #Cross-corr---------------------------------------------------------------- 
    [xcf,lid] = acu.xC_fft(ts,T);
    
    xc_p      = xcf[:,:,2:T-1];
    xc_p      = np.flip(xc_p,axis = 2); #positive-lag xcorrs
    xc_n      = xcf[:,:,T+1:-1];        #negative-lag xcorrs
    
    
    #Monster Equation---------------------------------------------------------
    wgt     = np.arange(nLg,1,-1)
    wgtm2   = np.tile((np.tile(wgt,[N,1])),[N,1]);
    wgtm3   = np.reshape(wgtm2,[N,N,np.size(wgt)]); #this is shit, eats all the memory!
    Tp      = T-1;
    
    # Da Equation!--------------------
     VarHatRho = (Tp*(1-rho**2)**2 ...
     +   rho**2  * sum(wgtm3 * (SumMat(ac**2,nLg)  +  xc_p**2 + xc_n**2),3)\         #1 2 4
     -   2*rho   * sum(wgtm3 * (SumMat(ac,nLg)     * (xc_p    + xc_n))  ,3)\         # 5 6 7 8
     +   2       * sum(wgtm3 * (ProdMat(ac,nLg)    + (xc_p   .* xc_n))  ,3))./(T**2);     
    
def SumMat(Y0):     
    T = np.shape(Y0)[0]
    N = np.shape(Y0)[1]
    
    F = (N*(N-1))/2
    SM0 = np.empty(N,N,T)
    
    for i in 
    
    
def ProdMat(Y0):    
    
def CorrMat(ts,Type='',**kwargs):
    """ Produce sample correlation matrices or Naively corrected z maps."""
    import numpy as np
    N   = np.shape(ts)[0];
    R = np.corrcoef(ts)
    R[range(N),range(N)] = 0;
    if Type=='' or Type=='rho':
        R[range(N),range(N)] = 0
        return(R)
    else:
        Z = np.arctanh(R)*np.sqrt(T-3)
        Z[range(N),range(N)] = 0
        return(Z)
    
 
################ TESTING PART:    
import scipy.io
#import numpy as np
V = '/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_124422_OnlyMTS.mat'
mat = scipy.io.loadmat(V)   
mts = mat['mts']
T = 1200
xDF_Calc(mts,T)