
clear
N=500;

nRlz=5000;

M=N;
B   = (eye(M)-ones(M)/M)./M;
ac_bnd=(sqrt(2)*erfinv(0.95))./sqrt(N);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
disp('=')
for i=1:nRlz
    c0    = corr(corrautocorr([0 0],0.8,0,N,1)');
    c(i) = c0(1,2);
end
mean(c)
v_null_MC_hc=var(atanh(c));
nf_v_null_MC_hc=var(c);

%ac1=[1 zeros(1,N-2)];
%ac2=[1 zeros(1,N-2)];

ts_null_hc = corrautocorr([0 0],0.8,0,N,1);
xAC=AC_fft(ts_null_hc,N);
%Threshold without sorting i.e. so, the far lags are ignored...
xAC(abs(xAC)<ac_bnd) = 0;
ac1 = xAC(1,:);
ac2 = xAC(2,:);

figure; hold on; imagesc(xAC)

%Threshold with sorting. 
%ac10 = xAC(1,1:end-1); ac10=ac10(abs(ac10)>ac_bnd); ac1=[ac10 zeros(1,N-(numel(ac10)+1))];
%ac20 = xAC(2,1:end-1); ac20=ac20(abs(ac20)>ac_bnd); ac2=[ac20 zeros(1,N-(numel(ac20)+1))];

v_null_C_hc=trace(toeplitz(ac1)*toeplitz(ac2))/(trace(toeplitz(ac1))*trace(toeplitz(ac2)));

v_null_D_hc=trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)));

[c0,v0]=HetBivCalc_fft(ts_null_hc,N);
sqrt(1./v0(1,2))

[v_null_D_hc v_null_C_hc v_null_MC_hc ]
sqrt(1./[v_null_D_hc v_null_C_hc v_null_MC_hc ])

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
disp('================================================')
clear c ac1 ac2 xAC
for i=1:nRlz
    c0    = corr(corrautocorr([0 0],0.2,0,N,1)');
    c(i) = c0(1,2);
end
mean(c)
v_null_MC_lc=var(atanh(c));
nf_v_null_MC_lc=var(c);

% ac1=[1 zeros(1,N-2)];
% ac2=[1 zeros(1,N-2)];

ts_null_lc = corrautocorr([0 0],0.2,0,N,1);
xAC=AC_fft(ts_null_lc,N);

%Threshold without sorting i.e. so, the far lags are ignored...
xAC(abs(xAC)<ac_bnd) = 0;
ac1 = xAC(1,:);
ac2 = xAC(2,:);

figure; hold on; imagesc(xAC)

%Threshold with sorting. 
%ac10 = xAC(1,1:end-1); ac10=ac10(abs(ac10)>ac_bnd); ac1=[ac10 zeros(1,N-(numel(ac10)+1))];
%ac20 = xAC(2,1:end-1); ac20=ac20(abs(ac20)>ac_bnd); ac2=[ac20 zeros(1,N-(numel(ac20)+1))];

v_null_C_lc=trace(toeplitz(ac1)*toeplitz(ac2))/(trace(toeplitz(ac1))*trace(toeplitz(ac2)));

v_null_D_lc=trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)));

[c0,v0]=HetBivCalc_fft(ts_null_lc,N);
sqrt(1./v0(1,2))

[v_null_D_lc v_null_C_lc v_null_MC_lc ]
sqrt(1./[v_null_D_lc v_null_C_lc v_null_MC_lc ])
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
disp('================================================')
clear c c ac1 ac2 xAC
for i=1:nRlz
    c0    = corr(corrautocorr([0 0],0.8,0.4,N,1)');
    c(i) = c0(1,2);
end
mean(c)
v_MC_AC1_hc=var(atanh(c));
nf_v_MC_AC1_hc=var(c);

% ac1 = [1 0.4 zeros(1,N-3)];
% ac2 = [1 0.4 zeros(1,N-3)];

ts_ac1_hc  = corrautocorr([0 0],0.8,0.4,N,1);
xAC=AC_fft(ts_ac1_hc,N);

%Threshold without sorting i.e. so, the far lags are ignored...
xAC(abs(xAC)<ac_bnd) = 0;
ac1 = xAC(1,:);
ac2 = xAC(2,:);

figure; hold on; imagesc(xAC)

%Threshold with sorting. 
%ac10 = xAC(1,1:end-1); ac10=ac10(abs(ac10)>ac_bnd); ac1=[ac10 zeros(1,N-(numel(ac10)+1))];
%ac20 = xAC(2,1:end-1); ac20=ac20(abs(ac20)>ac_bnd); ac2=[ac20 zeros(1,N-(numel(ac20)+1))];

v_C_AC1_hc=trace(toeplitz(ac1)*toeplitz(ac2))/(trace(toeplitz(ac1))*trace(toeplitz(ac2)));

v_D_AC1_hc=trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)));

[c0,v0]=HetBivCalc_fft(ts_ac1_hc,N);
sqrt(1./v0(1,2))

[v_D_AC1_hc v_C_AC1_hc v_MC_AC1_hc]
sqrt(1./[v_D_AC1_hc v_C_AC1_hc v_MC_AC1_hc])
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
disp('================================================')
clear c ac1 ac2 xAC 
for i=1:nRlz
    c0   = corr(corrautocorr([0 0],0.8,[0.6 .3],N,1)');
    c(i) = c0(1,2);
end
mean(c)
v_MC_AC3_hc=var(atanh(c));
nf_v_MC_AC3_hc=var(c);

% ac1=[1 [0.6 .4 .2] zeros(1,N-5)];
% ac2=[1 [0.6 .4 .2] zeros(1,N-5)];

ts_ac3_hc  = corrautocorr([0 0],0.8,[0.6 .3],N,1);
xAC=AC_fft(ts_ac3_hc,N);

%Threshold without sorting i.e. so, the far lags are ignored...
xAC(abs(xAC)<ac_bnd) = 0;
ac1 = xAC(1,:);
ac2 = xAC(2,:);

figure; hold on; imagesc(xAC)

%Threshold with sorting. 
%ac10 = xAC(1,1:end-1); ac10=ac10(abs(ac10)>ac_bnd); ac1=[ac10 zeros(1,N-(numel(ac10)+1))];
%ac20 = xAC(2,1:end-1); ac20=ac20(abs(ac20)>ac_bnd); ac2=[ac20 zeros(1,N-(numel(ac20)+1))];

v_C_AC3_hc=trace(toeplitz(ac1)*toeplitz(ac2))/(trace(toeplitz(ac1))*trace(toeplitz(ac2)));

v_D_AC3_hc=trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)));

[c0,v0]=HetBivCalc_fft(ts_ac3_hc,N);
sqrt(1./v0(1,2))

[v_D_AC3_hc v_C_AC3_hc v_MC_AC3_hc]
sqrt(1./[v_D_AC3_hc v_C_AC3_hc v_MC_AC3_hc])

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
disp('================================================')
clear c ac1 ac2 xAC 
for i=1:nRlz
    c0   = corr(corrautocorr([0 0],0.8,[-0.4 0.4],N,1)');
    c(i) = c0(1,2);
end
mean(c)
v_MC_AC3_hc=var(atanh(c));
nf_v_MC_AC3_hc=var(c);

% ac1=[1 [0.6 .4 .2] zeros(1,N-5)];
% ac2=[1 [0.6 .4 .2] zeros(1,N-5)];

ts_ac2neg_hc  = corrautocorr([0 0],0.8,[-0.4 0.4],N,1);
xAC=AC_fft(ts_ac2neg_hc,N);

%Threshold without sorting i.e. so, the far lags are ignored...
xAC(abs(xAC)<ac_bnd) = 0;
ac1 = xAC(1,:);
ac2 = xAC(2,:);

%figure; hold on; imagesc(xAC)

%Threshold with sorting. 
%ac10 = xAC(1,1:end-1); ac10=ac10(abs(ac10)>ac_bnd); ac1=[ac10 zeros(1,N-(numel(ac10)+1))];
%ac20 = xAC(2,1:end-1); ac20=ac20(abs(ac20)>ac_bnd); ac2=[ac20 zeros(1,N-(numel(ac20)+1))];

v_C_AC3_hc=trace(toeplitz(ac1)*toeplitz(ac2))/(trace(toeplitz(ac1))*trace(toeplitz(ac2)));

v_D_AC3_hc=trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)));

[c0,v0]=HetBivCalc_fft(ts_ac2neg_hc,N);
sqrt(1./v0(1,2))

[v_D_AC3_hc v_C_AC3_hc v_MC_AC3_hc]
sqrt(1./[v_D_AC3_hc v_C_AC3_hc v_MC_AC3_hc])
