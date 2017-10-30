clear

edfest = @(v) (1./v)+1;
c2t = @(C,edf) C.*sqrt((edf-2)./(1-C.^2));
bias= @(old,new) (old-new)./old*100;

disp('=================================')
addpath(genpath('/Users/sorooshafyouni/Home/GitClone/HetBiv'))

ii=1;
N    = 100;
%for ii=1:1000
    clearvars -except N c_null Curb t_Clif v_Dut v_Cliff b_Clif t_Dut b_Dut t_null ii bias c2t edfest
    ts  = corrautocorr([0 0],0.8,0,N,1);
    %ts  = ts-mean(ts,2);
    ts  = detrend(ts')';
    %ts   = randn(2,N);
    xAC = AC_fft(ts,N);
    ac1 = xAC(1,:); 
    ac2 = xAC(2,:);

    disp('CORR----------')
    cc=corr(ts');
    c_null(ii)=cc(1,2);
    c=cc(1,2)
    clearvars -except N c ts ac1 ac2 v_null Curb bias c2t v_Dut v_Cliff edfest tnull ii c_null t_Clif b_Clif t_Dut b_Dut t_null
    disp('TRUESHIT--------------------------')
    tnull=c2t(c,N)
    t_null(ii)=tnull;
    disp('Clliford--------------------------')
    v   = trace(toeplitz(ac1)*toeplitz(ac2))/(trace(toeplitz(ac1))*trace(toeplitz(ac2)));
    v_Cliff(ii)=v
    % B = eye(N);
    % v   = trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)))
    t   = corr(ts').*sqrt(edfest(v))./sqrt(1-corr(ts').^2);
    t_Clif(ii) = t(1,2)
    b_Clif(ii) = bias(tnull,t(1,2));
    
    disp('Dutilleul-------------------------')
    clearvars -except N ts ac1 ac2 v_null Curb bias c2t v_Dut v_Cliff edfest tnull ii c_null t_Clif b_Clif t_Dut b_Dut t_null
    B   = (eye(N)-ones(N)/N)./N;
    v   = trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)));
    v_Dut(ii)=v
    t   = corr(ts').*sqrt(edfest(v))./sqrt(1-corr(ts').^2);
    t_Dut(ii) = t(1,2)
    b_Dut(ii) = bias(tnull,t(1,2));

%     disp('=============')
%     disp('Naive Z--------------------------')
%     clearvars -except N ts ac1 ac2 bias c2t edfest tnull
%     %v   = HetBivCalc_fft(ts,N);
%     z   = atanh(corr(ts')).*sqrt(N);
%     turthz = z(1,2);
% 
%     %disp('HetBiv Z--------------------------')
%     clearvars -except N ts ac1 ac2 bias c2t edfest turthz
%     cf    = HetBivCalc_fft(ts,N,'method','CFx');
%     z     = atanh(corr(ts')).*sqrt(N/cf(1,2));
%     bz(i) = bias(turthz,z(1,2));
%end

% figure; 
% histogram(t_Clif,30,'Normalization','probability'); hold on
% histogram(t_Dut,30,'Normalization','probability'); hold on
% histogram(t_null,30,'Normalization','probability'); hold on