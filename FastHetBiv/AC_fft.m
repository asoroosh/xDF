function [xAC]=AC_fft(Y,L,varargin)
%[xAC]=AC_fft(Y,L,varargin)
% Super fast full-lag AC calculation of multi-dimention matrices. The
% function exploits fft to estimate the autocorrelations. 
%
%%%%INPUTS
%   Y:      A matrix of size IxT containing the I time series and T length.
%   L:      Time series length
%   
%   To get a double sided AC, add 'two-sided' as an input
%%%%OUTPUTS
%   xAC:    IxT-1 matrix of full-lag autocorrelations. If 'two-sided'
%           switched, then the matrix is Ix2(T-1). 
%_________________________________________________________________________
% Soroosh Afyouni, Ox, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if size(Y,2)~=L
    Y=Y';
end

Y=Y-mean(Y,2); 
%only works on >2016 Matlabs, but faster!
% if <2016, use Y=Y-repmat(mean(Y,2),1,L) instead.

nfft    = 2.^nextpow2(2*L-1);
yfft    = fft(Y,nfft,2);

xAC     = ifft(yfft.*conj(yfft),[],2)./sum(abs(Y).^2,2);
if sum(strcmpi(varargin,'two-sided'))
   xAC  = [xAC(:,end-L+2:end) ; xAC(:,1:L)];
else
    xAC=xAC(:,1:L);
end

