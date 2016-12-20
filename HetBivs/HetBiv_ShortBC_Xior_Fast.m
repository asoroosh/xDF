function [BCF]=HetBiv_ShortBC_Xior_Fast(ts1,ts2,ndpr,howfar)
% function [BCF]=HetBiv_ShortBC_Xior_Fast(ts1,ts2,ndpr,howfar)
% Thie function estimates a correction factor to calculate a Effective 
% Degree of Freedom (EDoF) for the case which the dependency between time
% series is limited to 0-lagged versions.
%
% NB! HetBiv_ShortBC_Xior and HetBiv_ShortBC_Xior_Fast calculate the same thing, but
% with remarkably different exec time.
% 
% Inputs:
% ts1 and ts2: timeseries of node 1 and node 2
% ndpr: number of time-points
% howfar:   curbing factor. Should be a value between 0 and 1:
%           howfar=round(1./TS) only takes the first lag
%           howfar=round((TS-1)./TS) consider all lags
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
a=autocorr(ts1,round(ndpr.*howfar));
b=autocorr(ts2,round(ndpr.*howfar));
BCF0=(ndpr+2.*sum(((ndpr-1):-1:1)'.*(a(2:end).*b(2:end))))./ndpr;

%% Second Part
% [ab,abl]=(crosscorr(ts1,ts2,round(ndpr.*howfar)));
% OneQCombCov=toeplitz(ab(abl>=0));
BCF1=(1./ndpr).*corr(a,b); %only one quarter of the dude!

%% And... the dude!
BCF=BCF0+BCF1;