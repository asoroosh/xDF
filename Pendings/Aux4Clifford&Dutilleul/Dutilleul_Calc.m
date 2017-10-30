function v=Dutilleul_Calc(ts,N)

M=N-1;
B      = (eye(M)-ones(M)/M)./M;
%ac_bnd = (sqrt(2)*erfinv(0.95))./sqrt(N);

xAC=AC_fft(ts,N);
ac1 = xAC(1,1:end-1);
ac2 = xAC(2,1:end-1); 

%v=trace(toeplitz(ac1)*toeplitz(ac2))/(trace(toeplitz(ac1))*trace(toeplitz(ac2)));

v=trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)));