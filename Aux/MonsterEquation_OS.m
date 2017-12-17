function [ASAt,CnR] = MonsterEquation_OS(rho,ACx,ACy,T)
% CALCULATE the variance for given rho ACx ACy and T
% ONLY use this for oracle simulations. 
%
% Soroosh Afyouni, Ox, 2017

Sigma_x  = MakeMeCovMat(ACx,T);
Sigma_y  = MakeMeCovMat(ACy,T);

Kx = sqrtm(Sigma_x);
Ky = sqrtm(Sigma_y);

as=max(diag(Kx*Ky'));

rho=rho*as;

Sigma_xy = (Kx*Ky').*rho/as; 

figure; imagesc(Sigma_xy)

%Monster Eq.
SigX2    = (rho.^2./2).* trace(Sigma_x ^2);
SigY2    = (rho.^2./2).* trace(Sigma_y ^2);
SigXSigY = trace(Sigma_x * Sigma_y);

CnR = SigXSigY/T^2;

SigXSigXY= -2.*rho     .* trace(Sigma_x * Sigma_xy);
SigYSigXY= -2.*rho     .* trace(Sigma_y * Sigma_xy);
SigXY2   = rho.^2.* trace(Sigma_xy^2)+trace(Sigma_xy^2);

ASAt=(SigX2+SigY2+SigXSigY+SigXY2+SigXSigXY+SigYSigXY)/T.^2;

%(1-rho.^2).^2/T


