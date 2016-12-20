function [BCF,BCF0,BCF1]=HetBiv_ShortBC_Fast(ts1,ts2,ndpr,howfar)
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
a=autocorr(ts1,ndpr-1);
b=autocorr(ts2,ndpr-1);
BCF0=(ndpr+2.*sum(((ndpr-1):-1:1)'.*(a(2:end).*b(2:end))))./ndpr;
%% Second Part
ns_ndpr=round(ndpr.*howfar);
ab=crosscorr(ts1,ts2,ns_ndpr-1);
BCF1=sum((ns_ndpr-abs((ns_ndpr-1):-1:(-1*ns_ndpr+1))').*(ab.^2))./ndpr;
%% And... the dude!
BCF=BCF0+BCF1;