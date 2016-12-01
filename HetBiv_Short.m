function [BCF]=HetBiv_Short(ts1,ts2,ndpr,howfar,figflag)
% function [BCF]=HetBiv_Short(ts1,ts2,ndpr,howfar)
% Thie function estimates a correction factor to calculate a Effective 
% Degree of Freedom (EDoF) for the case of *independence* between two time
% series.
%
% NB! HetBiv_Short and HetBiv_Short_Fast calculate the same thing, but
% with remarkably different exec time.
% 
% Inputs:
% ts1 and ts2: timeseries of node 1 and node 2
% ndpr: number of time-points
% howfar:   curbing factor. Should be a value between 0 and 1:
%           howfar=round(1./TS) only takes the first lag
%           howfar=round((TS-1)./TS) consider all lags
% figflag:  should be set to 1 if you want to see toeplitz and cov
%           matrices, otherwise 0.
%
% Output:
% BCF:  Correction Factor. i.e. EDoF=DoF/BCF.
%
%

if size(ts1,1)~=(ndpr)  
    disp(['TS1 timeseries transposed!'])
    ts1=ts1';
end
if size(ts2,1)~=(ndpr)
    disp(['TS2 timeseries transposed!'])
    ts2=ts2';
end

[a,al]=autocorr(ts1,round(ndpr.*howfar));
ta=toeplitz(a);
[b,bl]=autocorr(ts2,round(ndpr.*howfar));
tb=toeplitz(b);
CombCov=ta*tb;
BCF=(1./ndpr).*trace(CombCov);

if figflag
    disp(['The dude is: ' num2str(BCF)])
    figure; hold on; 
    subplot(1,3,1); hold on; 
    imagesc(al,al,ta);colorbar; axis image
    subplot(1,3,2); hold on; 
    imagesc(bl,bl,tb);colorbar; axis image    
    subplot(1,3,3); hold on;
    imagesc(al,al,CombCov);colorbar; axis image
end