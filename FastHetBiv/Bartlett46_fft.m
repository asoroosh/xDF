function [BCF,BCFA]=Bartlett46_fft(Y,L)
% [BCF,BCFA]=Bartlett46_fft(Y,L)
% Super-fats calculation of Bartlett correction factors. 
%
%%%INPUTS:
%   Y: IxT or TxI BOLD voxel-wise/pacellated time series 
%   L: Length of the time series! Just to check the dimensions are all
%   sorted!
%%%OUTPUTS:
%   BCF:  In the classic 1935 Bartlett's Correction Factor. \hat{N}=N/BCF
%   BCFA: Is the Richardson & Clifford's modification for \rho>0. \hat{N}=N/BCFA
%
%%%DEPENDENCIES:
%   AC_fft.m 
% 
%%%REFERENCES:
%   Bretherton et al, Journal of Climate, 1999, p2004
%   Robert Haining, Geographical Analysis, 1991
%   Richardson & Clifford, 1991, p300
%_________________________________________________________________________
% Soroosh Afyouni, NISOx.org, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________
%
if size(Y,2)~=L
    Y=Y'; %IxT
end

xAC = AC_fft(Y,L);
xAC = xAC(:,2); % AC lag-1

xAC=xAC*xAC';

BCF=(1+xAC)./(1-xAC);
BCFA=(1+xAC)./((1-xAC).*(1-corr(Y').^2).^2);