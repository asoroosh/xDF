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


##############################################################################
########################## START OF xDF_Calc #################################
##############################################################################

def xDF_Calc(ts,T,\
             method      = 'truncate',\
             methodparam = 'adaptive',\
             verbose     = True,\
             TV          = True,\
             copy        = True):
    
    # -------------------------------------------------------------------------
##### READ AND CHECK 0---------------------------------------------------------
    
    #if not verbose: blockPrint()
    
    if copy: #Make sure you are not messing around with the original time series
        ts = ts.copy()
    
    if np.shape(ts)[1]!=T:
        if verbose: print('xDF::: Input should be in IxT form, the matrix was transposed.')
        ts = np.transpose(ts)
    
    N  = np.shape(ts)[0];
    
    ts_std = np.std(ts,axis=1,ddof=1)
    ts     = ts/np.transpose(np.tile(ts_std,(T,1)));   #standardise
    print('xDF_Calc::: Time series standardised by their standard deviations.')
    # READ AND CHECK 0---------------------------------------------------------    
    # -------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------
##### Estimate xC and AC ------------------------------------------------------
    
    #Corr----------------------------------------------------------------------
    rho,znaive    = CorrMat(ts,T)
    rho           = np.round(rho,7)
    znaive        = np.round(znaive,7)
    #Autocorr------------------------------------------------------------------
    [ac,CI] = AC_fft(ts,T); 
    ac   = ac[:,1:T-1]; #The last element of ACF is rubbish, the first one is 1, so why bother?!
    nLg  = T-2;  
    
    #Cross-corr---------------------------------------------------------------- 
    [xcf,lid] = xC_fft(ts,T);
    
    xc_p      = xcf[:,:,1:T-1];
    xc_p      = np.flip(xc_p,axis = 2); #positive-lag xcorrs
    xc_n      = xcf[:,:,T:-1];        #negative-lag xcorrs

    # -------------------------------------------------------------------------    
##### Start of Regularisation--------------------------------------------------
    if method.lower()=='tukey':
        if methodparam=='':
            M = np.sqrt(T)
        else: M = methodparam
        if verbose: print('xDF_Calc::: AC Regularisation: Tukey tapering of M = ' + str(int(np.round(M))))
        ac   = tukeytaperme(ac,nLg,M)
        xc_p = tukeytaperme(xc_p,nLg,M)
        xc_n = tukeytaperme(xc_n,nLg,M)
        
        #print(np.round(ac[0,0:50],4))
        
    elif method.lower()=='truncate':
        if type(methodparam)==str:    #Adaptive Truncation
            if methodparam.lower()!='adaptive':
                raise ValueError('What?! Choose adaptive as the option, or pass an integer for truncation')
            if verbose: print('xDF_Calc::: AC Regularisation: Adaptive Truncation')         
            [ac,bp] = shrinkme(ac,nLg)
            #truncate the cross-correlations, by the breaking point found from the ACF. (choose the largest of two)
            for i in np.arange(N):
                for j in np.arange(N):
                    maxBP        = np.max([bp[i],bp[j]])
                    xc_p[i,j,:]  = curbtaperme(xc_p[i,j,:],nLg,maxBP,verbose=False)
                    xc_n[i,j,:]  = curbtaperme(xc_n[i,j,:],nLg,maxBP,verbose=False)
        elif type(methodparam) == int: #Npne-Adaptive Truncation
            if verbose: print('xDF_Calc::: AC Regularisation: Non-adaptive Truncation on M = ' + str(methodparam))         
            ac    = curbtaperme(ac,nLg,methodparam)
            xc_p  = curbtaperme(xc_p,nLg,methodparam)
            xc_n  = curbtaperme(xc_n,nLg,methodparam)
            
        else: raise ValueError('xDF_Calc::: methodparam for truncation method should be either str or int.')
    # Start of Regularisation--------------------------------------------------
    # -------------------------------------------------------------------------
    
    
    # -------------------------------------------------------------------------
##### Start of the Monster Equation--------------------------------------------
    # -------------------------------------------------------------------------
    
    wgt     = np.arange(nLg,0,-1)
    wgtm2   = np.tile((np.tile(wgt,[N,1])),[N,1]);
    wgtm3   = np.reshape(wgtm2,[N,N,np.size(wgt)]); #this is shit, eats all the memory!
    Tp      = T-1;
    
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

    # -----------------------------------------------------------------------
    # End of the Monster Equation--------------------------------------------
    # -----------------------------------------------------------------------
    
    # -----------------------------------------------------------------------
##### Truncate to Theoritical Variance --------------------------------------
    TV_val = (1-rho**2)**2/T;
    TV_val[range(N),range(N)] = 0
        
    idx_ex = np.where(VarHatRho < TV_val)
    NumTVEx = (np.shape(idx_ex)[1])/2;
    #print(NumTVEx)
    
    if NumTVEx>0 and TV:
        if verbose: print('Variance truncation is ON.')
        # Assuming that the variance can *only* get larger in presence of autocorrelation.  
        VarHatRho[idx_ex] = TV_val[idx_ex];
        #print(N)
        #print(np.shape(idx_ex)[1])
        FGE = N*(N-1)/2
        if verbose: print('xDF_Calc::: ' + str(NumTVEx) + ' (' + str(round((NumTVEx/FGE)*100,3)) + '%) edges had variance smaller than the textbook variance!')
    else: 
        if verbose: print('xDF_Calc::: NO truncation to the theoritical variance.')
    # Sanity Check:
    #        for ii in np.arange(NumTVEx):
    #            print( str( idx_ex[0][ii]+1 ) + '  ' + str( idx_ex[1][ii]+1 ) )


    # -------------------------------------------------------------------------        
#####Start of Statistical Inference -------------------------------------------
    
    # Well, these are all Matlab and pretty useless -- copy pasted them just in case though...
    #Pearson's turf -- We don't really wanna go there, eh?
    #rz      = rho./sqrt((ASAt));     %abs(ASAt), because it is possible to get negative ASAt!
    #r_pval  = 2 * normcdf(-abs(rz)); %both tails
    #r_pval(1:nn+1:end) = 0;          %NaN screws up everything, so get rid of the diag, but becareful here. 
    
    #Our turf--------------------------------    
    rf      = np.arctanh(rho)
    sf      = VarHatRho/((1-rho**2)**2)    #delta method; make sure the N is correct! So they cancel out.
    rzf     = rf/np.sqrt(sf)
    f_pval  = 2 * sp.norm.cdf(-abs(rzf))  #both tails
    
     # diagonal is rubbish;
    VarHatRho[range(N),range(N)] = 0
    f_pval[range(N),range(N)]    = 0         #NaN screws up everything, so get rid of the diag, but becareful here. 
    rzf[range(N),range(N)]       = 0
    
    #End of Statistical Inference ---------------------------------------------
    # -------------------------------------------------------------------------    
    
    #if not verbose: enablePrint()
    
    xDFOut = {'p':f_pval,\
              'z':rzf,\
              'znaive': znaive,\
              'v':VarHatRho,\
              'TV':TV_val,\
              'TVExIdx':idx_ex}
    
    return xDFOut

##############################################################################
########################## END OF xDF_Calc ###################################
##############################################################################
    
# Disable verbose
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore verbose
def enablePrint():
    sys.stdout = sys.__stdout__    
    
    
