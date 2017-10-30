function v=Clifford_Calc_Trimmed(ts,N)

M=N-1;
%B      = (eye(M)-ones(M)/M)./M;
ac_bnd = (sqrt(2)*erfinv(0.95))./sqrt(N);

xAC=AC_fft(ts,N);
ac100 = xAC(1,1:end-1); ac10=ac100(abs(ac100)>ac_bnd); ac1=[ac10 zeros(1,N-(numel(ac10)+1))];
ac200 = xAC(2,1:end-1); ac20=ac200(abs(ac200)>ac_bnd); ac2=[ac20 zeros(1,N-(numel(ac20)+1))];

v=trace(toeplitz(ac1)*toeplitz(ac2))/(trace(toeplitz(ac1))*trace(toeplitz(ac2)));
%v=trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)));