function [ASAt,CnR] = MonsterEquation_OS(rho,ACx,ACy,T)
% CALCULATE the variance for given rho ACx ACy and T
% ONLY use this for oracle simulations. 

Sigma_x  = MakeMeCovMat(ACx,T);
Sigma_y  = MakeMeCovMat(ACy,T);

Kx = sqrtm(Sigma_x);
Ky = sqrtm(Sigma_y);

d   = diag(Kx*Ky');

as  = min(d);
rho = rho.*as;

Sigma_xy = rho.*((Kx*Ky')./as);

% figure; plot(d); hold on; plot(diag((Kx*Ky')./as)); hold off; 
% figure; 
% subplot(3,1,1); hold on; axis tight; 
% imagesc(Sigma_x)
% subplot(3,1,2); hold on; axis tight; 
% imagesc(Sigma_y)
% subplot(3,1,3); hold on; axis tight; 
% imagesc(Sigma_xy)
%Monster Eq -----------------------------------------------
SigX2    =  trace(Sigma_x ^2);
SigY2    =  trace(Sigma_y ^2);
SigXSigY =  trace(Sigma_x * Sigma_y);

SigXY2    = trace(Sigma_xy^2);

SigXSigXY = trace(Sigma_x * Sigma_xy);
SigYSigXY = trace(Sigma_y * Sigma_xy);

ASAt=((rho.^2./2).* SigX2... 
    + (rho.^2./2) .* SigY2...
    + rho.^2      .* SigXY2...
    - 2.*rho      .* SigXSigXY...
    - 2.*rho      .* SigYSigXY... 
    + SigXSigY...
    + SigXY2)/T.^2;
%----------------------------------------------------------
CnR = SigXSigY/T^2;
%(1-rho.^2).^2/T


