function t=mvnr_CholDecomp(mu,sig,ndp)
%t=corrautocorr(mu,sig,ndp)
%
%Multivariate Normal Distribution via Cholesky Decomposition
%
% Similar to MATLAB built-in mvnrnd.m 
%
% SA, Ox, 2017

z = randn(numel(mu),ndp); %Z
t = chol(sig)'*z+mu';      %X=M+BZ for X~N(mu,BB')

