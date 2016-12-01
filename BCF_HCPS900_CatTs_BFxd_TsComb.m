%[BCF,info]=BCF_HCPS900_CatTs_BFxd_TsComb(ts1,ts2,ndpr,NumStd,HowFarTs)
%Estimates a correction factor for sample size of timeseries affected by
%autocorrelation
%
%ts1,ts2:   two input timeseries - if you want to calc a correction factor for
%           only one timeseries, set both the same. 
%ndpr:      #datapoints
%NumStd:    Number of std for adjusting the boundry in 'smart' case
%HowFarTs:  A value between 0 to 1 detemining how long of a timeseries
%           should be considered in 'dumb' case. 
%
%
%   Consideration: 
%   1)  This is function is to estimate a correction factor according to 
%       unbiased variance of *Autocorrelation Function*.
%       In other words, it is mathematically correct version of what 
%       Fox et al 2005 and Van Dijk 2012 *wrongly* suggest for correcting 
%       the DoF. Their proposal in basically wrong because:
%           A) The Bartlett's 1945 fails to poropose a unbiased variance estimation
%              of correlation coefficients and BCF=\sum_{k=-\infty}^{+\infty}\rho_{k} is for
%              variance of ACF not correlation coefficients. 
%           B) Even if it was for correlation coefficients, in practice we don't
%              have -infty of ACF, therefore the pratical form should be as follow:
%
%              BCF=1+(2\sum_{k=1}^{K-1}\rho_{k})
%
%   2) This faction is for *local* estimation of correction factor, whereas,
%      Fox et al 2005 and Van Dijk 2012, *globally* estimate this
%      correction factor. For global estimation of correction factor.
%
%
%SA-2016-UoW

function [BCF,info]=BCF_HCPS900_CatTs_BFxd_TsComb(ts1,ts2,ndpr,NumStd,HowFarTs)

% NumStd=2;

if size(ts1,1)~=(ndpr)  
    disp(['TS1 timeseries transposed!'])
    ts1=ts1';
%     size(ts1)
end
if size(ts2,1)~=(ndpr)
    disp(['TS2 timeseries transposed!'])
    ts2=ts2';
%     size(ts2)
end

if nargin==5
    ExceededP=[];smartflag=0;ExcdIdx=[];bounds2=[];bounds1=[];
    
    ACFtmp1=(autocorr(ts1',round(ndpr.*HowFarTs)));
    ACFtmp2=(autocorr(ts2',round(ndpr.*HowFarTs)));

    bcf_s=1+(2*sum(ACFtmp1(2:end).*ACFtmp2(2:end)));
    clear ACtmp2 ACtmp1
elseif nargin==4
%     disp('smart')
    smartflag=1;
    [ACFtmp1,~,bounds1]=autocorr(ts1',ndpr-1,[],NumStd); %bound was calculated by the 2 numSTD as defualt
    [ACFtmp2,~,bounds2]=autocorr(ts2',ndpr-1,[],NumStd); %bound was calculated by the 2 numSTD as defualt
    
    ACFtmp1=ACFtmp1(2:end); ACFtmp2=ACFtmp2(2:end);
%     if bounds1(1)~=bounds2(1); warning('ACF bounds are different!'); end; 
    ExcdIdx=sqrt(abs(ACFtmp1).*abs(ACFtmp2))>bounds1(1);
    ExceededP=sum(ExcdIdx)./ndpr;
 
    bcf_s=1+(2.*sum((ACFtmp1(ExcdIdx)).*(ACFtmp2(ExcdIdx))));
    
    clear SigACF ACF bounds ACFtmp1 ACFtmp2
else
    error('Something is wrong with inputs!')
end

info.ExceededP=ExceededP;
info.Smart=smartflag;
info.ExcdIdx=ExcdIdx;
info.bounds1=bounds1;
info.bounds2=bounds2;

% info=[];

BCF=bcf_s;
    