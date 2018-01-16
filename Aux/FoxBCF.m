function [Z,BCF] = FoxBCF(Y,T)
% The way that Fox and VanDjk estimate the BCF; i.e. Integral [sci] of ACFs 
% Fox et al 2005 and VanDijek 200X
%
% SA, Ox, 2018

if size(Y,1)~=T
    Y=Y'; %IxT
end

ac        = AC_fft(Y,T);
ac(:,end) = [];
BCF       = sum(ac.^2,2);
BCF       = mean(BCF);
Z         = atanh(corr(Y)).*sqrt(T./BCF-3);