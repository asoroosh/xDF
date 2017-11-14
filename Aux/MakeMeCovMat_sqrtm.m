function [sigCmat]=MakeMeCovMat_sqrtm(sigC,ndp)

if size(sigC,1)==numel(sigC); sigC=sigC'; end;
db0     = [fliplr(sigC) 1 sigC];
sigC0   = spdiags(ones(ndp,1)*db0,(-numel(sigC):numel(sigC)),ndp,ndp);
sigCmat = full(sigC0);

