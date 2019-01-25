#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 16:41:29 2019


Many of these functions were copy pasted from bctpy package:
    https://github.com/aestrivex/bctpy
under GNU V3.0:
    https://github.com/aestrivex/bctpy/blob/master/LICENSE
    


@author: sorooshafyouni
University of Oxford, 2019
"""

import numpy as np
import statsmodels.stats.multitest as smmt
import scipy.stats as sp
import matplotlib.pyplot as plt


#def sub2ind(DaShape, X, Y):
#    """ DaShape: [#rows #columns]
#        returns idx of given X & Y
#        SA, Ox, 2018 """
#    return X*DaShape[1] + Y

def OLSRes(YOrig,RG,T,copy=True):
    """ 
    Or how to deconfound stuff!
    For regressing out stuff from your time series, quickly and nicely!
    SA,Ox,2019
    """
    if copy:
        YOrig = YOrig.copy()
        
    if np.shape(YOrig)[0]!=T or np.shape(RG)[0]!=T:
        raise ValueError('The Y and the X should be TxI format.')
        
    #demean anyways!
    mRG = np.mean(RG,axis=0)
    RG = RG-np.tile(mRG,(T,1));
    #B       = np.linalg.solve(RG,YOrig) # more stable than pinv
    invRG   = np.linalg.pinv(RG)    
    B       = np.dot(invRG,YOrig)
    Yhat    = np.dot(RG,B) # find the \hat{Y}
    Ydeconf = YOrig-Yhat #get the residuals -- i.e. cleaned time series
    
    return Ydeconf

def issymmetric(W):
    """Check whether a matrix is symmetric"""
    return((W.transpose() == W).all())
    
def SumMat(Y0,T,copy=True): 
    """
    Parameters
    ----------
    Y0 : a 2D matrix of size TxN

    Returns
    -------
    SM : 3D matrix, obtained from element-wise summation of each row with other 
         rows.    
         
    SA, Ox, 2019     
    """    
    
    if copy:
        Y0 = Y0.copy()    
    
    if np.shape(Y0)[0]!=T:
        print('SumMat::: Input should be in TxN form, the matrix was transposed.')
        Y0 = np.transpose(Y0)
    
    N = np.shape(Y0)[1]
    Idx = np.triu_indices(N)
    #F = (N*(N-1))/2
    SM = np.empty([N,N,T])
    for i in np.arange(0,np.size(Idx[0])-1):
        xx = Idx[0][i]
        yy = Idx[1][i]
        SM[xx,yy,:] = (Y0[:,xx]+Y0[:,yy]);
        SM[yy,xx,:] = (Y0[:,yy]+Y0[:,xx]);
    
    return SM

def ProdMat(Y0,T,copy=True):
    """
    Parameters
    ----------
    Y0 : a 2D matrix of size TxN

    Returns
    -------
    SM : 3D matrix, obtained from element-wise multiplication of each row with 
         other rows. 
         
    SA, Ox, 2019       
    """    
    
    if copy:
        Y0 = Y0.copy()

    if np.shape(Y0)[0]!=T:
        print('ProdMat::: Input should be in TxN form, the matrix was transposed.')
        Y0 = np.transpose(Y0)

    N = np.shape(Y0)[1]
    Idx = np.triu_indices(N)
    #F = (N*(N-1))/2
    SM = np.empty([N,N,T])
    for i in np.arange(0,np.size(Idx[0])-1):
        xx = Idx[0][i]
        yy = Idx[1][i]
        SM[xx,yy,:] = (Y0[:,xx]*Y0[:,yy]);
        SM[yy,xx,:] = (Y0[:,yy]*Y0[:,xx]);
    
    return SM

def CorrMat(ts,T,method='rho',copy=True):
    """ 
    Produce sample correlation matrices 
    or Naively corrected z maps.
    """

    if copy:
        ts = ts.copy()

    if np.shape(ts)[1]!=T:
        print('xDF::: Input should be in IxT form, the matrix was transposed.')
        ts = np.transpose(ts)

    N = np.shape(ts)[0];
    R = np.corrcoef(ts)
    
    Z = np.arctanh(R)*np.sqrt(T-3)
    
    R[range(N),range(N)] = 0
    Z[range(N),range(N)] = 0
    
    return R, Z
    
    #R[range(N),range(N)] = 0
#    if method.lower()=='rho':
#        R[range(N),range(N)] = 0
#        return(R)
#    elif method.lower()=='naive':
#        Z = np.arctanh(R)*np.sqrt(T-3)
#        Z[range(N),range(N)] = 0
#        return(Z)
#    else: print('Choose between either Type==`Naive` or Type==`rho`')


def stat_threshold(Z,mce='fdr_bh',a_level=0.05,side='two',copy=True):
    """
    Threshold z maps
    
    Parameters
    ----------
    
    mce: multiple comparison error correction method, should be
    among of the options below. [defualt: 'fdr_bh']. 
    The options are from statsmodels packages:
        
        `b`, `bonferroni` : one-step correction
        `s`, `sidak` : one-step correction
        `hs`, `holm-sidak` : step down method using Sidak adjustments
        `h`, `holm` : step-down method using Bonferroni adjustments
        `sh`, `simes-hochberg` : step-up method  (independent)
        `hommel` : closed method based on Simes tests (non-negative)
        `fdr_i`, `fdr_bh` : Benjamini/Hochberg  (non-negative)
        `fdr_n`, `fdr_by` : Benjamini/Yekutieli (negative)
        'fdr_tsbh' : two stage fdr correction (Benjamini/Hochberg)
        'fdr_tsbky' : two stage fdr correction (Benjamini/Krieger/Yekutieli)
        'fdr_gbs' : adaptive step-down fdr correction (Gavrilov, Benjamini, Sarkar)
    """
    
    if copy:
        Z = Z.copy()
    
    if side=='one':
        sideflag = 1
    elif side=='two' or 'double':
        sideflag = 2
    
    Idx = np.triu_indices(Z.shape[0],1)
    Zv = Z[Idx]
    
    Pv = sp.norm.cdf(-np.abs(Zv))*sideflag
            
    [Hv,adjpvalsv] = smmt.multipletests(Pv,method = mce)[:2]    
    adj_pvals = np.zeros(Z.shape)
    Zt = np.zeros(Z.shape)
            
    Zv[np.invert(Hv)] = 0 
    Zt[Idx] = Zv
    Zt = Zt + Zt.T; 
    
    adj_pvals[Idx] = adjpvalsv   
    adj_pvals = adj_pvals + adj_pvals.T; 
    
    adj_pvals[range(Z.shape[0]),range(Z.shape[0])] = 0
                
    return Zt, binarize(Zt), adj_pvals 

def RemoveNeg(Mats,copy=True):
    """ quickly remove negative values"""
    if copy:
        Mats = Mats.copy()
    Mats[Mats<0] = 0
    return(Mats)
    

class MatManParamError(RuntimeError):
    pass


def threshold_absolute(W, thr, copy=True):
    '''
    This function thresholds the connectivity matrix by absolute weight
    magnitude. All weights below the given threshold, and all weights
    on the main diagonal (self-self connections) are set to 0.

    If copy is not set, this function will *modify W in place.*

    Parameters
    ----------
    W : np.ndarray
        weighted connectivity matrix
    thr : float
        absolute weight threshold
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : np.ndarray
        thresholded connectivity matrix
    '''
    if copy:
        W = W.copy()
    np.fill_diagonal(W, 0)  # clear diagonal
    W[W < thr] = 0  # apply threshold
    return W


def threshold_proportional(W, p, copy=True):
    '''
    This function "thresholds" the connectivity matrix by preserving a
    proportion p (0<p<1) of the strongest weights. All other weights, and
    all weights on the main diagonal (self-self connections) are set to 0.

    If copy is not set, this function will *modify W in place.*

    Parameters
    ----------
    W : np.ndarray
        weighted connectivity matrix
    p : float
        proportional weight threshold (0<p<1)
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : np.ndarray
        thresholded connectivity matrix

    Notes
    -----
    The proportion of elements set to 0 is a fraction of all elements
    in the matrix, whether or not they are already 0. That is, this function
    has the following behavior:

    >> x = np.random.random((10,10))
    >> x_25 = threshold_proportional(x, .25)
    >> np.size(np.where(x_25)) #note this double counts each nonzero element
    46
    >> x_125 = threshold_proportional(x, .125)
    >> np.size(np.where(x_125))
    22
    >> x_test = threshold_proportional(x_25, .5)
    >> np.size(np.where(x_test))
    46

    That is, the 50% thresholding of x_25 does nothing because >=50% of the
    elements in x_25 are aleady <=0. This behavior is the same as in BCT. Be
    careful with matrices that are both signed and sparse.
    '''
    from .miscellaneous_utilities import teachers_round as round

    if p > 1 or p < 0:
        raise MatManParamError('Threshold must be in range [0,1]')
    if copy:
        W = W.copy()
    n = len(W)						# number of nodes
    np.fill_diagonal(W, 0)			# clear diagonal

    if np.allclose(W, W.T):				# if symmetric matrix
        W[np.tril_indices(n)] = 0		# ensure symmetry is preserved
        ud = 2						# halve number of removed links
    else:
        ud = 1

    ind = np.where(W)					# find all links

    I = np.argsort(W[ind])[::-1]		# sort indices by magnitude

    en = int(round((n * n - n) * p / ud))		# number of links to be preserved

    W[(ind[0][I][en:], ind[1][I][en:])] = 0  # apply threshold
    #W[np.ix_(ind[0][I][en:], ind[1][I][en:])]=0

    if ud == 2:						# if symmetric matrix
        W[:, :] = W + W.T						# reconstruct symmetry

    return W


def weight_conversion(W, wcm, copy=True):
    '''
    W_bin = weight_conversion(W, 'binarize');
    W_nrm = weight_conversion(W, 'normalize');
    L = weight_conversion(W, 'lengths');

    This function may either binarize an input weighted connection matrix,
    normalize an input weighted connection matrix or convert an input
    weighted connection matrix to a weighted connection-length matrix.

    Binarization converts all present connection weights to 1.

    Normalization scales all weight magnitudes to the range [0,1] and
    should be done prior to computing some weighted measures, such as the
    weighted clustering coefficient.

    Conversion of connection weights to connection lengths is needed
    prior to computation of weighted distance-based measures, such as
    distance and betweenness centrality. In a weighted connection network,
    higher weights are naturally interpreted as shorter lengths. The
    connection-lengths matrix here is defined as the inverse of the
    connection-weights matrix.

    If copy is not set, this function will *modify W in place.*

    Parameters
    ----------
    W : NxN np.ndarray
        weighted connectivity matrix
    wcm : str
        weight conversion command.
        'binarize' : binarize weights
        'normalize' : normalize weights
        'lengths' : convert weights to lengths (invert matrix)
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : NxN np.ndarray
        connectivity matrix with specified changes

    Notes
    -----
    This function is included for compatibility with BCT. But there are
    other functions binarize(), normalize() and invert() which are simpler to
    call directly.
    '''
    if wcm == 'binarize':
        return binarize(W, copy)
    elif wcm == 'normalize':
        return normalize(W, copy)
    elif wcm == 'lengths':
        return invert(W, copy)
    else:
        raise NotImplementedError('Unknown weight conversion command.')


def binarize(W, copy=True):
    '''
    Binarizes an input weighted connection matrix.  If copy is not set, this
    function will *modify W in place.*

    Parameters
    ----------
    W : NxN np.ndarray
        weighted connectivity matrix
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : NxN np.ndarray
        binary connectivity matrix
    '''
    if copy:
        W = W.copy()
    W[W != 0] = 1
    return W


def normalize(W, copy=True):
    '''
    Normalizes an input weighted connection matrix.  If copy is not set, this
    function will *modify W in place.*

    Parameters
    ----------
    W : np.ndarray
        weighted connectivity matrix
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : np.ndarray
        normalized connectivity matrix
    '''
    if copy:
        W = W.copy()
    W /= np.max(np.abs(W))
    return W


def invert(W, copy=True):
    '''
    Inverts elementwise the weights in an input connection matrix.
    In other words, change the from the matrix of internode strengths to the
    matrix of internode distances.

    If copy is not set, this function will *modify W in place.*

    Parameters
    ----------
    W : np.ndarray
        weighted connectivity matrix
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : np.ndarray
        inverted connectivity matrix
    '''
    if copy:
        W = W.copy()
    E = np.where(W)
    W[E] = 1. / W[E]
    return W


def autofix(W, copy=True):
    '''
    Fix a bunch of common problems. More specifically, remove Inf and NaN,
    ensure exact binariness and symmetry (i.e. remove floating point
    instability), and zero diagonal.


    Parameters
    ----------
    W : np.ndarray
        weighted connectivity matrix
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : np.ndarray
        connectivity matrix with fixes applied
    '''
    if copy:
        W = W.copy()

    # zero diagonal
    np.fill_diagonal(W, 0)

    # remove np.inf and np.nan
    W[np.logical_or(np.where(np.isinf(W)), np.where(np.isnan(W)))] = 0

    # ensure exact binarity
    u = np.unique(W)
    if np.all(np.logical_or(np.abs(u) < 1e-8, np.abs(u - 1) < 1e-8)):
        W = np.around(W, decimal=5)

    # ensure exact symmetry
    if np.allclose(W, W.T):
        W = np.around(W, decimals=5)

    return W

def density_und(CIJ):
    '''
    Density is the fraction of present connections to possible connections.

    Parameters
    ----------
    CIJ : NxN np.ndarray
        undirected (weighted/binary) connection matrix

    Returns
    -------
    kden : float
        density
    N : int
        number of vertices
    k : int
        number of edges

    Notes
    -----
    Assumes CIJ is undirected and has no self-connections.
            Weight information is discarded.
    '''
    n = len(CIJ)
    k = np.size(np.where(np.triu(CIJ).flatten()))
    kden = k / ((n * n - n) / 2)
    return kden, n, k