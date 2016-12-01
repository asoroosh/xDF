%[BCF35,BCF35A]=BCF35nArbab_HCPS900_CatTs_TsComb(ts1,ts2,r,ndpr)
%
%   This function *locally* (i.e. per edge) estimates the correction factor. 
%
%Inputs:
%   ts1 and ts2:    Time series from node 1 and node 2
%   r:  correlation coefficient between ts1 and ts2
%   ndpr:   number of data-points in ts1 and ts2
%           ndpr of node1 and node2 should be equal
%
%Outputs:
%   BCF35:  AR(1) Bartlett's 'correction factor' (Bartlett 1935)
%   BCF35A: Unbiased version of BCF (Arbabshirani et al 2015)
%
%
% SA-2016-UoW
function [BCF35,BCF35A]=BCF35nArbab_HCPS900_CatTs_TsComb(ts1,ts2,r,ndpr)

% ndpr=1200;

if size(ts1,1)~=(ndpr)  
    disp(['TS1 timeseries transposed!'])
    ts1=ts1';
end
if size(ts2,1)~=(ndpr)
    disp(['TS2 timeseries transposed!'])
    ts2=ts2';
end

% ts1=nets_demean(ts1);
ts1=ts1-mean(ts1);
ARonetmp_0=sum(ts1(1:end-1).*ts1(2:end))/sum(ts1.*ts1);
ts2=ts2-mean(ts2);
ARonetmp_1=sum(ts2(1:end-1).*ts2(2:end))/sum(ts2.*ts2);

BCF35=((1+(ARonetmp_0.*ARonetmp_1))./(1-(ARonetmp_0.*ARonetmp_1)));
BCF35A=((1-r^2)^2).*(1+(ARonetmp_0.*ARonetmp_1))./(1-(ARonetmp_0.*ARonetmp_1));

