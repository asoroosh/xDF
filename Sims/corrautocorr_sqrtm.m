function [t]=corrautocorr_sqrtm(mu,sigR,sigC,ndp)
% t = corrautocorr(mu,sig,ndp)
%   Draws random numbers from Matrix Normal Distribution
%%%INPUTS
%   mu  : Expected values. A vector of size 1xP where P is number of time 
%         series required
%   sigR: Covariance matrix (between time series; i.e. Correlations).
%         If only two univariate time series required a scalar (indicating 
%         the off-diagonal elements) suffice. 
%   SigC: Covariance of columns (within time series; i.e. Serial-correlations)
%         By default the code assumes you need a similar AC structure for
%         all time-series. If you need different structures for each time 
%         series, feed SigC with a 3D matrix of NxNxP where N is number of 
%         data-points i.e. time series p follows AC structur on p^th layer 
%         of the SigC matrix. This takes much longer than univeral AC 
%         structure though!
%   ndp : Number of data-points (N)
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

if isvector(sigC)
   if size(sigC,1) == numel(sigC); sigC = sigC'; end;
   db0   = [fliplr(sigC) 1 sigC];
   sigC0 = spdiags(ones(ndp,1)*db0,(-numel(sigC):numel(sigC)),ndp,ndp);
   sigC  = full(sigC0);
end
if isscalar(sigR) %assumes it is a lazy way of saying only two time series!
    assert(numel(mu) == 2, 'sigR should be a square matrix!')
    sigR0 = ones(2)*sigR;
    sigR0(1:2+1:end) = 1; %variance of x and y is 1.
    sigR  = sigR0;
end

numt    = numel(mu);
z       = randn(numt,ndp); % X~N(0,1)
[CsigR] = sqrtm(sigR)'; %the transpose is ineffective, however it matches the Chol version.
MVDisk  = CsigR*z;

if size(sigC,3)>1
    assert(size(sigC,3) == numt,'Number of layers in SigC should match the number of time series.')
    for l = 1:size(sigC,3)
        sigC_tmp  = sigC(:,:,l);
        CsigC_tmp = sqrtm(sigC_tmp)';
        t(l,:)    = mu(l)+MVDisk(l,:)*CsigC_tmp;
        clear *_tmp
    end
else
    CsigC = sqrtm(sigC);
    t     = mu'+MVDisk*CsigC';    % Y=M+AXB for Y~MN(mu,AA',BB')
    t     = real(t);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Old version 
% z = randn(numel(mu),ndp);            
%t_test = sqrtm(sigR)'*z*sqrtm(sigC)'+mu';    % Y=M+AXB for Y~MN(mu,AA',BB')
%corr(t_test')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
