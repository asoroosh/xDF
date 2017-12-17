function [ASA,CnR]=MonsterEquation(ts,T)

if size(ts,1) ~= 2 || size(ts,2) ~= T
    error('PFF!')
end

%ts  = dtrend(ts);
ts  = ts./std(ts,[],2); %standardise
%Corr----------------------------------------------------------------------
c   = corr(ts');
rho = c(1,2);
%Autocorr------------------------------------------------------------------
[ac] = AC_fft(ts,T);
%Tukey Tappering----------------------------------------
%    Mac = T/5; %round(sqrt(T));
%    ac(1,:) = tukeytaperme(ac(1,:),Mac,1);
%    ac(2,:) = tukeytaperme(ac(2,:),Mac,1);
%Shrinkage------------------------------------------------
ac(1,:) = shrinkme(ac(1,:));
ac(2,:) = shrinkme(ac(2,:));
%Curbing------------------------------------------------
% Mac = 9;%T/5; %round(sqrt(T));
% ac  = curbtaperme(ac,Mac);
%------------------------------------------------
ac_x = ac(1,1:T-1);
ac_y = ac(2,1:T-1);

Sigma_x  = toeplitz(ac_x);
Sigma_y  = toeplitz(ac_y);

%figure; imagesc(Sigma_x)
%Cross-corr---------------------------------------------------------------- 
 [xcf,lags]  = crosscorr(ts(1,:),ts(2,:),T-1);

 ut_CC = fliplr(xcf(2:T)); %rhos sit hear, so later, we use -1 diag.
 lt_CC = xcf(T:end-1);
 
%Tukey Tapering the CROSS correlations----------------------------   
%    Mcc = T/5; %;round(sqrt(T));
%    ut_CC = tukeytaperme(ut_CC,Mcc,rho);
%    lt_CC = tukeytaperme(lt_CC,Mcc,rho);
%Shrinkage------------------------------------------------
ut_CC = shrinkme(ut_CC);
lt_CC = shrinkme(lt_CC);

plot(ut_CC)
%Curbing---------------------------------------------------------
% Mcc   = 9; %T/5; %round(sqrt(T));
% ut_CC = curbtaperme(ut_CC,Mcc);
% lt_CC = curbtaperme(lt_CC,Mcc);
%Smoothing-------------------------------------------------------
 %KrnlSz = 2;
 %Sigma_xy = conv2(Sigma_xy, ones(KrnlSz)/KrnlSz^2,'same');
 %Sigma_xy = eye(T-1).*rho; % no cross-correlation; only 0lag corr
 %figure; imagesc(Sigma_xy) 
%--------------------------------------------------------------------------
Sigma_xy = (triu(toeplitz(ut_CC))+tril(toeplitz(lt_CC),-1))'; 
figure; imagesc(Sigma_xy)
% figure; 
% subplot(3,2,1); hold on;
% title('Sigma_{xy}')
% imagesc(Sigma_xy,[-1 1])
% disp(['Sigma_{xy}: ' num2str(trace(Sigma_xy))])
% 
% subplot(3,2,2); hold on; 
% title('Sigma_x x Sigma_y')
% imagesc(Sigma_x * Sigma_y,[-1 1])
% disp(['Sigma_{x}Sigma_{y}: ' num2str(trace(Sigma_x * Sigma_y)) ]);
% 
% subplot(3,2,3); hold on; 
% title('Sigma_x x Sigma_{xy}')
% imagesc(Sigma_x * Sigma_xy,[-1 1])
% disp(['Sigma_{x}Sigma_{xy}: ' num2str(trace(Sigma_x * Sigma_xy)) ]);
% 
% subplot(3,2,4); hold on; 
% title('Sigma_y x Sigma_{xy}')
% imagesc(Sigma_y * Sigma_xy,[-1 1])
% disp(['Sigma_{y}Sigma_{xy}: ' num2str(trace(Sigma_y * Sigma_xy))])
% 
% subplot(3,2,5); hold on; 
% title('Sigma x')
% imagesc(Sigma_x,[-1 1])
% disp(['Sigma_{x}: ' num2str(trace(Sigma_x))])
% 
% subplot(3,2,6); hold on; 
% title('Sigma y')
% imagesc(Sigma_y,[-1 1])
% disp(['Sigma_{y}: ' num2str(trace(Sigma_y))])

CnRe = trace(Sigma_x * Sigma_y); % Clifford and Richardson

ASA=((rho.^2./2).* trace(Sigma_x ^2        )...
    +rho.^2     .* trace(Sigma_xy^2        )...
    -2.*rho     .* trace(Sigma_x * Sigma_xy)...
    +(rho.^2./2).* trace(Sigma_y ^2        )...
    -2.*rho     .* trace(Sigma_y * Sigma_xy)...
    + CnRe...
    + trace(Sigma_xy^2))./T^2;          
%--------------------------------------------------------------------------          
CnR=CnRe./T^2; %just to check how bad the others are doing!


function ct_ts=curbtaperme(ts,M)
  M          = round(M);
  msk        = zeros(size(ts));
  msk(:,1:M) = 1;
  ct_ts      = msk.*ts;
  
function srnkd_ts=shrinkme(ts)
L = numel(ts);
bnd = (sqrt(2)*erfinv(0.95))./sqrt(L);
idx = find(abs(ts)>bnd);
isit       = abs(ts)>bnd & (1:L);
where2stop = find(isit==0);
where2stop = where2stop(1);
srnkd_ts   = curbtaperme(ts,where2stop);

function tt_ts=tukeytaperme(ts,M,intv)
  if ~exist('intv','var'); intv = 1; warning('Oi!'); end;
  M          = round(M);
  tt_ts      = zeros(size(ts));
  tt_ts(:,1) = intv;
  tt_ts(2:M) = (1+cos([2:M].*pi./M))./2.*ts(2:M);
