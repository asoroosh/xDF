function [t]=corrautocorr(mu,sigR,sigC,ndp,verboseflag)
% t = corrautocorr(mu,sig,ndp)
%   Draws from Matrix Normal Distribution by Cholesky Decomposition
%
%   NB! If any of 2nd/3rd input is not positive semidefinite matrix, then it
%   automatically converts them to a nearest PSD matrix 
%   following Highman, 1988, Linear Algebra and its Applications.
%
%%%INPUTS
%   mu  : Expected values. A vector of size 1xP where P is number of time 
%         series required
%
%   sigR: Covariance matrix (between time series; i.e. Correlations).
%         If only two univariate time series required a scalar (indicating 
%         the off-diagonal elements) suffice. 
%
%   SigC: Covariance of columns (within time series; i.e. Serial-correlations)
%         By default the code assumes you need a similar AC structure for
%         all time-series. If you need different structures for each time 
%         series, feed SigC with a 3D matrix of NxNxP where N is number of 
%         data-points i.e. time series p follows AC structur on p^th layer 
%         of the SigC matrix. This takes much longer than univeral AC 
%         structure though!
%
%   ndp : Number of data-points (N)
%
%   verboseflag [optional] : if you want to see the warning for non-PSD 
%                            input matrix [default: 1] 
%
%%%NOTE
%   If you need to use a non-PSD matrix for very large simulations use
%   function "nearestSPD" outside of "corrautocorr" for once. 
%
%%%OUTPUTS
%   t   : Simulated time series
%
%%%EXAMPLE
%   This produces two autocorrelated time series of 1000 length whith 0.6 
%   correlation between them 
%
%   ts=corrautocorr([0 0],0.6,[0.6 0.3 0.2],1000);
%   corr(ts')
%   ac1=autocorr(ts(1,:)), ac2=autocorr(ts(2,:)) 
%
%%%REFERENCES:
%   Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%   under Serial Correlation
%   Soroosh Afyouni, Stephen M. Smith & Thomas E. Nichols
%   2018
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________


%Exp sin damped function to model the AC behaviour -- Should be checked
%with the real data in case it is too off!
%   Amp=1; f=10; atcon=7; phi=0; N=1200;
%   t=0:round(1/(N/4),5):1;
%   f=f*2*pi;
%   y=Amp*cos(f*t + phi).*exp(-atcon*t); plot(y)

%ts=corrautocorr([0 0],0.5,[0.8 0.4 0.5 0.3 0.3],N); xAC=AC_fft(ts,N); sac=xAC./sqrt(sum(xAC(:,1:N/4).^2,2)/N); p=fdr_bh(2.*normcdf(-abs(sac)))

if ~exist('verboseflag','var'); verboseflag=1; end;
if isvector(sigC)
   if size(sigC,1)==numel(sigC); sigC=sigC'; end;
   db0   = [fliplr(sigC) 1 sigC];
   sigC0 = spdiags(ones(ndp,1)*db0,(-numel(sigC):numel(sigC)),ndp,ndp);
   sigC  = full(sigC0);
end
if isscalar(sigR) %assumes it is a lazy way of saying only two time series!
    assert(numel(mu)==2, 'sigR should be a square matrix!')
    sigR0 = ones(2)*sigR;
    sigR0(1:2+1:end) = 1; %variance of x and y is 1.
    sigR = sigR0;
end
numt  = numel(mu);
z     = randn(numt,ndp); % X~N(0,1)


[CsigR,psdflag] = chol(sigR); %sigR=(chol(abs(sigR)).*sign(sigR))'; %that is stupid!
if psdflag
    if verboseflag; warning('Found the nearest PSD matrix for Corr.'); end;
    sigR  = nearestSPD(sigR); %in case they are not SPD
    CsigR = chol(sigR);
end

MVDisk = CsigR'*z;

if size(sigC,3)>1
    assert(size(sigC,3)==numt,'Number of layers in SigC should match the number of time series.')
    for l=1:size(sigC,3)
        sigC_tmp  = sigC(:,:,l);
        sigC_tmp  = nearestSPD(sigC_tmp); 
        CsigC_tmp = chol(sigC_tmp)';
        t(l,:)    = mu(l)+MVDisk(l,:)*CsigC_tmp;
        clear *_tmp
    end
else
    %disp('Similar AC structure across time series.')
    %imagesc(sigC)
    %CsigC = (chol(abs(sigC)).*sign(sigC))'; %that is stupid!
    [CsigC,psdflag]=chol(sigC);
    if psdflag
        clear CsigC
        sigC  = nearestSPD(sigC); 
        %figure; imagesc(sigC);
        if verboseflag; warning('Found the nearest PSD matrix for AC.'); end;
        CsigC = chol(sigC);
    end
    t     = mu'+MVDisk*CsigC';    % Y=M+AXB for Y~MN(mu,AA',BB')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Old version 
% z = randn(numel(mu),ndp);            
%t_test = sqrtm(sigR)'*z*sqrtm(sigC)'+mu';    % Y=M+AXB for Y~MN(mu,AA',BB')
%corr(t_test')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%EXTERNAL FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ahat = nearestSPD(A)
% From Higham: "The nearest symmetric positive semidefinite matrix in the
% Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
% where H is the symmetric polar factor of B=(A + A')/2."
%
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%-------------------------------------------------------------------------
%Copyright (c) 2013, John D'Errico 
% All rights reserved.
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% * Redistributions of source code must retain the above copyright 
% notice, this list of conditions and the following disclaimer. 
% * Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in 
% the documentation and/or other materials provided with the distribution

if nargin ~= 1
  error('Exactly one argument must be provided.')
end
[r,c] = size(A);
if r ~= c
  error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
  Ahat = eps;
  return
end
B = (A + A')/2;
[~,Sigma,V] = svd(B);
H    = V*Sigma*V';
Ahat = (B+H)/2;
Ahat = (Ahat + Ahat')/2;
p    = 1;
k    = 0;
while p ~= 0
  [~,p]  = chol(Ahat);
  k = k + 1;
  if p  ~= 0
    mineig = min(eig(Ahat));
    Ahat   = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
  end
end
