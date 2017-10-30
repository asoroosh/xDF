function [v,vv]=Clifford_Calc(ts,N)

%M=N-1;
%B      = (eye(M)-ones(M)/M)./M;
%ac_bnd = (sqrt(2)*erfinv(0.95))./sqrt(N);

[xAC,~,xACov]=AC_fft(ts,N);

ac1   = xAC(1,1:end-1);
ac2   = xAC(2,1:end-1);

acov1 = xACov(1,1:end-1);
acov2 = xACov(2,1:end-1);

v=trace(toeplitz(acov1)*toeplitz(acov2))/(trace(toeplitz(acov1))*trace(toeplitz(acov2)));

vv=trace(toeplitz(ac1)*toeplitz(ac2))/(trace(toeplitz(ac1))*trace(toeplitz(ac2)));
%v=trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)));