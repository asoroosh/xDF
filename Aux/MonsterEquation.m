function [SAS,CR]=MonsterEquation(ts,T)

if size(ts,1) ~= 2 || size(ts,2) ~= T
    error('PFF!')
end

ts  = ts./std(ts,[],2); %standardise
%Corr
c   = corr(ts');
rho = c(1,2)

%Curbing Factors-----------------------------------------------------------
% M    = round(sqrt(T));
% N    = M; 

%M = T-1;
M = round(sqrt(T));
N = T-1;

%Autocorr------------------------------------------------------------------
ac   = AC_fft(ts,T);

ac_tmp    = zeros(size(ac));
ac_tmp(:,1) = 1;
ac_tmp(:,2:M) = (1+cos([2:M].*pi./M))./2.*ac(:,2:M);

ac = ac_tmp; clear *_tmp;
figure; plot(ac')
ac_x = ac(1,1:T-1);
ac_y = ac(2,1:T-1);

Sigma_x  = toeplitz(ac_x);
Sigma_y  = toeplitz(ac_y);
figure; imagesc(Sigma_x)
%Cross-corr----------------------------------------------------------------
% xcf   = crosscorr(ts(1,:),ts(2,:),M);
% 
% ut_CC = fliplr(xcf(1:M+1)) .*[ones(1,N) zeros(1,M+1-N)];
% lt_CC = xcf(M+1:end)       .*[ones(1,N) zeros(1,M+1-N)];
% 
% Sigma_xy = triu(toeplitz(ut_CC))+tril(toeplitz(lt_CC),-1);
%--------------------------------------------------------------------------

Sigma_xy = eye(T-1).*rho; % no cross-correlation; only 0lag corr
figure; imagesc(Sigma_xy)
CRe = trace(Sigma_x * Sigma_y); %Clifford and Richardson

SAS=((rho^2/2)* trace(Sigma_x^2         )...
    +rho^2    * trace(Sigma_xy^2        )...
    -2*rho    * trace(Sigma_x * Sigma_xy)...
    +(rho^2/2)* trace(Sigma_y^2         )...
    -2*rho    * trace(Sigma_y * Sigma_xy)...
              + CRe                      ...
              + trace(Sigma_xy^2        ))./T^2;          
%--------------------------------------------------------------------------          
CR=CRe./T^2;


% SAS=((rho^2/2)* trace(Sigma_x^2         ) ./v_x.^2   ...
%     +rho^2    * trace(Sigma_xy^2        ) ./(v_x*v_y)...
%     -2*rho    * trace(Sigma_x * Sigma_xy) ./(v_x*s_x*s_y)  ...
%     +(rho^2/2)* trace(Sigma_y^2         ) ./v_y.^2  ...
%     -2*rho    * trace(Sigma_y * Sigma_xy) ./(v_y*s_y*s_x)  ...
%               + CRe                       ./(v_x*v_y)...
%               + trace(Sigma_xy^2        ) ./(v_x*v_y))./T^2;