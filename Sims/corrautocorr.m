function t=corrautocorr(mu,sigR,sigC,ndp)
% t = corrautocorr(mu,sig,ndp)
% Multivariate Normal Distribution via Cholesky Decomposition
%
%%%INPUTS
%   mu  : Expected values
%   sigR: Covariance of rows (between time series; i.e. Correlations)
%   SigC: Covariance of columns (within time series; i.e. Serial-correlations)
%   ndp : Number of data-points
%
%%%OUTPUTS
%   t   : Simulated time series
%
% srafyouni@gmail.com
% SA, Ox, 2017

if isvector(sigC)
   db0   = [sigC 1 sigC];
   sigC0 = spdiags(ones(ndp,1)*db0,(-numel(sigC):numel(sigC)),ndp,ndp);
   sigC  = full(sigC0);
end

if isscalar(sigR) %assumes it is a lazy way of saying only two time series!
    assert(numel(mu)==2, 'sigR should be a square matrix!')
    sigR0 = ones(2)*sigR;
    sigR0(1:2+1:end) = 1; %variance of x and y is 1.
    sigR = sigR0;
end

z = randn(numel(mu),ndp);           % X~MN(mu,I,I)
t = chol(sigR)'*z*chol(sigC)'+mu';    % Y=M+AXB for Y~MN(mu,AA',BB')
