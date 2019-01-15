#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 12:25:21 2019

@author: sorooshafyouni
University of Oxford, 2019
"""

import scipy.io
from xDF import shrinkme
import numpy as np

V = '/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_124422_OnlyMTS.mat'
mat = scipy.io.loadmat(V)   
mts = mat['mts']
T = 1200

AC,bnd = AC_fft(mts,T)

print(bnd)

shrnk_ts,bp = shrinkme(AC[:2,:],T)