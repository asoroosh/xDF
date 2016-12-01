% function [BCF]=HetBiv_Short_Fast(ts1,ts2,ndpr,howfar)
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
%
% Output:
% BCF:  Correction Factor. i.e. EDoF=DoF/BCF.

function [BCF]=HetBiv_Short_Fast(ts1,ts2,ndpr,howfar)
if size(ts1,1)~=(ndpr)  
    disp(['TS1 timeseries transposed!'])
    ts1=ts1';
end
if size(ts2,1)~=(ndpr)
    disp(['TS2 timeseries transposed!'])
    ts2=ts2';
end

a=autocorr(ts1,round(ndpr.*howfar));
b=autocorr(ts2,round(ndpr.*howfar));
BCF=(ndpr+2.*sum(((ndpr-1):-1:1).*(a(2:end).*b(2:end))))./ndpr;