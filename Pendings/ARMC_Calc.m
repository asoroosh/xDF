function v=ARMC_Calc(ts,N,ACn)

if ~exist('ACn','var'); ACn=3; end;  

xAC = AC_fft(ts,N);
xAC = mean(xAC);
xAC = xAC(1:ACn);

C=corr(ts'); C=C(1,2);

nRlz=5000;
for i=1:nRlz
    c0    = corr(corrautocorr([0 0],C,xAC,N,1)');
    c(i) = c0(1,2);
end

v=var(atanh(c));