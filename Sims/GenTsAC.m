
function ts=GenTsAC(ar1,ndpr)
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

ts(1)=randn(1); 
for t=2:ndpr 
    ts(t)=ts(t-1)*ar1+randn(1); 
end