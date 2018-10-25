function [sigCmat]=MakeMeCovMat_sqrtm(sigC,ndp)
% this should be removed later! It was only for testing porpusing!
%
%function ts=GenTsAC(ar1,ndpr)
% Generates very dumb non-white noises 
% ar1:  The first-order autoregressive coefficient
% ndpr: Number of time points
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
if size(sigC,1)==numel(sigC); sigC=sigC'; end;
db0     = [fliplr(sigC) 1 sigC];
sigC0   = spdiags(ones(ndp,1)*db0,(-numel(sigC):numel(sigC)),ndp,ndp);
sigCmat = full(sigC0);

