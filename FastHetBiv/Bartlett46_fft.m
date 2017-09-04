function [BCF,BCFA]=Bartlett46_fft(Y,L)

if size(Y,2)~=L
    Y=Y'; %IxT
end
I = size(Y,1);

xAC = AC_fft(Y,L);
xAC = xAC(:,2); % AC lag-1

xAC=xAC*xAC';

BCF=(1-xAC)./(1+xAC);
BCFA=(1-xAC)./((1+xAC).*(1-corr(Y').^2));