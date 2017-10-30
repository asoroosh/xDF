function [sigCmat,CsigC]=MakeMeCovMat(sigC,ndp,verboseflag)

if ~exist('verboseflag','var'); verboseflag=1; end;

if size(sigC,1)==numel(sigC); sigC=sigC'; end;
db0   = [fliplr(sigC) 1 sigC];
sigC0 = spdiags(ones(ndp,1)*db0,(-numel(sigC):numel(sigC)),ndp,ndp);
sigCmat  = full(sigC0);

[CsigC,psdflag]=chol(sigCmat);
if psdflag
    clear CsigC
    if verboseflag; warning('Found the nearest PSD matrix for AC.'); end;
    sigCmat  = nearestSPD(sigCmat); %in case they are not SPD
    %figure; imagesc(sigC)
    CsigC = chol(sigCmat);
end

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