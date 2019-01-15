#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 11:56:37 2019

@author: sorooshafyouni
University of Oxford, 2019
"""

import scipy.io
from xDF import *
import numpy as np

V = '/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_124422_OnlyMTS.mat'
mat = scipy.io.loadmat(V)   
mts = mat['mts']
T = 1200

print('+++++++ xDF without regularisation::: +++++++++++++++++++++++++++++++')
xDFOut_TVOff = xDF_Calc(mts,T,method='',TV=False)
xDFOut_TVOn  = xDF_Calc(mts,T,method='',TV=True)

print('+++++++ xDF without truncation::: +++++++++++++++++++++++++++++++++++')
xDFOut_tna   = xDF_Calc(mts,T,method='truncate',methodparam = int(T/4),verbose=True)

print('+++++++ xDF with ADAPTIVE truncation::: ++++++++++++++++++++++++++++++')
xDFOut_ta    = xDF_Calc(mts,T,method='truncate',methodparam = 'adaptive',verbose=True)

Z = stat_threshold(xDFOut_ta['z'],mce='b')[0]
IDX = np.where(Z==0)
print(np.size(IDX[0])/2)

print('+++++++ xDF without tapering::: +++++++++++++++++++++++++++++++++++++')
xDFOut_tt    = xDF_Calc(mts,T,method='tukey',methodparam = np.sqrt(T),verbose=True)

#tt = tukeytaperme(ac,T,np.sqrt(T))
#[ac,CI] = AC_fft(mts,T)
#[sh, bp] = shrinkme(ac[1:T],T-1)