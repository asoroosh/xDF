function [ASA,CnR]=MonsterEquation(ts,T)

if size(ts,1) ~= 2 || size(ts,2) ~= T
    error('PFF!')
end

ts  = ts./std(ts,[],2); %standardise
%Corr
c   = corr(ts');
rho = c(1,2);

%Curbing Factors-----------------------------------------------------------
% M    = round(sqrt(T));
% N    = M; 

N = T-1;
%Autocorr------------------------------------------------------------------
ac   = AC_fft(ts,T);

%Tukey Tappering----------------------------------------
% M             = round(sqrt(T));
% ac_tmp        = zeros(size(ac));
% ac_tmp(:,1)   = 1;
% ac_tmp(:,2:M) = (1+cos([2:M].*pi./M))./2.*ac(:,2:M);
% ac            = ac_tmp; clear *_tmp;
% 
% plot(ac')

%Curbing------------------------------------------------
%ac_tmp  = ac;
%M       = round(T/5);
%acpvals = zeros(size(ac_tmp));
%acpvals(:,1:M) = 1;
%acpvals = ones(size(ac_tmp));
%ac      = acpvals.*ac_tmp;
%clear ac_tmp;

%figure; plot(ac')

ac_x = ac(1,1:T-1);
ac_y = ac(2,1:T-1);

Sigma_x  = toeplitz(ac_x);
Sigma_y  = toeplitz(ac_y);

%figure; imagesc(Sigma_x)
%Cross-corr---------------------------------------------------------------- 
 xcf   = crosscorr(ts(1,:),ts(2,:),T-1);

 % These two lines are for curbed versions:
 % ut_CC = fliplr(xcf(1:L)) .*[ones(1,G) zeros(1,M+1-G)];
 % lt_CC = xcf(M+1:end)       .*[ones(1,G) zeros(1,M+1-G)];
 
 ut_CC = fliplr(xcf(2:T)); %rhos sit hear, so later, we use -1 diag.
 lt_CC = xcf(T:end-1);
 
 %Tukey Tapering the CROSS correlations-------------------------
 %plot(ut_CC); hold on; plot(lt_CC); 
 
%  M = round(9*sqrt(T));
%  
%  ut_CC_tmp      = zeros(size(ut_CC));
%  ut_CC_tmp(:,1) = 1;
%  ut_CC_tmp(2:M) = (1+cos([2:M].*pi./M))./2.*ut_CC(2:M);
%  
%  lt_CC_tmp      = zeros(size(lt_CC));
%  lt_CC_tmp(:,1) = 1;
%  lt_CC_tmp(2:M) = (1+cos([2:M].*pi./M))./2.*lt_CC(2:M);
%  
%  ut_CC = ut_CC_tmp;
%  lt_CC = lt_CC_tmp;
 
 %plot(ut_CC,'linewidth',1.3); hold on; plot(lt_CC,'linewidth',1.3); legend({'2:T','T:End','TT 2:T','TT T:End'})
 
 %Smoothing-------------------------------------------------------
 %KrnlSz = 2;
 %Sigma_xy = conv2(Sigma_xy, ones(KrnlSz)/KrnlSz^2,'same');
 %Sigma_xy = eye(T-1).*rho; % no cross-correlation; only 0lag corr
 %figure; imagesc(Sigma_xy) 
%--------------------------------------------------------------------------

Sigma_xy = triu(toeplitz(ut_CC))+tril(toeplitz(lt_CC),-1); 

% figure; imagesc(Sigma_xy)

CnRe = trace(Sigma_x * Sigma_y); % Clifford and Richardson

% disp('----')
% disp('-2.*rho.*trace(Sigma_x * Sigma_xy)')
% -2.*rho     .* trace(Sigma_x * Sigma_xy)
% trace(Sigma_x * Sigma_xy)
% disp('-2.*rho.*trace(Sigma_y * Sigma_xy)')
% -2.*rho     .* trace(Sigma_y * Sigma_xy)
% trace(Sigma_y * Sigma_xy)
% disp('----')

ASA=((rho.^2./2).* trace(Sigma_x^2         )...
    +rho.^2     .* trace(Sigma_xy^2        )...
    -2.*rho     .* trace(Sigma_x * Sigma_xy)...
    +(rho.^2./2).* trace(Sigma_y^2         )...
    -2.*rho     .* trace(Sigma_y * Sigma_xy)...
                 + CnRe                     ...
                 + trace(Sigma_xy^2        ))./T^2;          
%--------------------------------------------------------------------------          
CnR=CnRe./T^2;


% SAS=((rho^2/2)* trace(Sigma_x^2         ) ./v_x.^2   ...
%     +rho^2    * trace(Sigma_xy^2        ) ./(v_x*v_y)...
%     -2*rho    * trace(Sigma_x * Sigma_xy) ./(v_x*s_x*s_y)  ...
%     +(rho^2/2)* trace(Sigma_y^2         ) ./v_y.^2  ...
%     -2*rho    * trace(Sigma_y * Sigma_xy) ./(v_y*s_y*s_x)  ...
%               + CRe                       ./(v_x*v_y)...
%               + trace(Sigma_xy^2        ) ./(v_x*v_y))./T^2;
















