function [v]=MCVar_Calc(x,y,L,nRlz)
% [v]=MCVar_Calc(x,y,L)
% Returns the MC variance of simulated time series
%
% SA, Ox, 2017

%assert(size(x)==size(y),'x and y should be the same size.')

if ~exist('nRlz','var'); nRlz=1000; end;

if size(x,1)~=L; x=x'; end
if size(y,1)~=L; y=y'; end

xAC=AC_fft([x,y],L);
xAC=xAC(:,2:5);

mxAC=mean(xAC); %there should be a better way!! Yeah, use different autocorrelation structures for each time series to corrautocorr.m

z=[];
for i=1:nRlz
    t=corrautocorr([0 0],corr(x,y),mxAC,L);
    z=[z;atanh(corr(t(1,:)',t(2,:)'))];
end
v=var(z);