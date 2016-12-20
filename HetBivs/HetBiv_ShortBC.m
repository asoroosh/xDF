function [BCF]=HetBiv_ShortBC(ts1,ts2,ndpr,howfar,figflag)
% function [BCF]=HetBiv_ShortBC_Fast(ts1,ts2,ndpr,howfar)
% Thie function estimates a correction factor to calculate a Effective 
% Degree of Freedom (EDoF) for the case of *dependency* between lagged
% versions of timeseries.
%
% NB! HetBiv_ShortBC and HetBiv_ShortBC_Fast calculate the same thing, but
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
%% First Part
[a,al]=autocorr(ts1,ndpr-1);
ta=toeplitz(a);
[b,bl]=autocorr(ts2,ndpr-1);
tb=toeplitz(b);
CombVar=ta*tb;
BCF0=(1./ndpr).*trace(CombVar);

%% Second Part
[ab,abl]=crosscorr(ts1,ts2,round(ndpr.*howfar));

OneQCombCov=toeplitz(ab(abl>=0));
BCF_correction=(1./ndpr).*trace(OneQCombCov^2); %only one quarter of the dude!

%% And... the dude!
BCF=BCF0+BCF_correction;

%% 
if figflag
    disp(['The dude is: ' num2str(BCF)])
    figure; hold on; 
    subplot(1,3,1); hold on; 
    imagesc(al,al,ta);colorbar; axis image
    subplot(1,3,2); hold on; 
    imagesc(bl,bl,tb);colorbar; axis image    
    subplot(1,3,3); hold on;
    imagesc(al,al,CombVar);colorbar; axis image
    
    figure; hold on; 
    subplot(1,3,1); hold on;
    plot(abl,ab)   
    subplot(1,3,2); hold on;
    imagesc(abl,abl,toeplitz(ab));axis image;colorbar 
    subplot(1,3,3); hold on;
    imagesc(abl(abl>=0),abl(abl>=0),OneQCombCov^2,[-.3,.3]);axis image;colorbar
end