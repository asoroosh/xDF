function [BCF,BCFA]=Bartlett46_fft(Y,L)
%
%
%%%REFERENCES:
%   BRETHERTON et al, Journal of Climate, 1999, p2004
%   Arbabshirani et al, Neuroimage, 2016, p 
%   Robert Haining, Geographical Analysis, 1991, p
%   Richardson & Clifford, 1991, p300
%
if size(Y,2)~=L
    Y=Y'; %IxT
end
%I = size(Y,1);

xAC = AC_fft(Y,L);
xAC = xAC(:,2); % AC lag-1

xAC=xAC*xAC';

BCF=(1+xAC)./(1-xAC);
BCFA=(1+xAC)./((1-xAC).*(1-corr(Y').^2).^2);