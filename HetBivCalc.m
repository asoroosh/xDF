function [CF,EDF,Stat]=HetBivCalc(X,Y,varargin)
% [CF]=HetBivCalc(x,y,varargin)
%   Calculates the effective degrees of freedom of correlation coefficient
%   of two time series X and Y. 
%
%%%INPUTS
%
%
%
%%%OUTPUTS
%
%
%
%%%EXAMPLES
%
%
%%%REFERENCE
%
%   Soroosh Afyouni, UoW, 2017
%   srafyouni@gmail.com
%
%
%


%Y1=randn(1000,1);
%Y2=randn(1000,1);
%Yt1=toeplitz(autocorr(Y1,length(Y1)-1)); 
%Yt2=toeplitz(autocorr(Y2,length(Y1)-1)); 
%trace(Yt1*Yt2)./length(Y1), 
%sum(sum(Yt1.*Yt2,2))./length(Y1)

ndpr = length(X);
%% Read

%
verbose=1; rho=[];
if sum(strcmpi(varargin,'verbose'))
   verbose = varargin{find(strcmpi(varargin,'verbose'))+1}; 
end

if sum(strcmpi(varargin,'Method'))
   CFm = varargin{find(strcmpi(varargin,'Method'))+1};
    if strcmpi(CFm,'BCF')
        CFmethod=1;
        bnd_cc=NaN;
    elseif strcmpi(CFm,'CF0')
        CFmethod=2;
        bnd_cc=NaN;
    elseif strcmpi(CFm,'CFx')
        CFmethod=3;
        if sum(strcmpi(varargin,'CCHowfar'))
           bnd_cc = varargin{find(strcmpi(varargin,'CCHowfar'))+1};
        else
           bnd_cc = 10; 
        end
    else
        error('Unknown BCF methods. Choose between BCF, BC0, BCx');
    end
else
   CFm='BCF';
   CFmethod=1;
   bnd_cc=NaN;
end
%
if sum(strcmpi(varargin,'Howfar'))
   bnd = varargin{find(strcmpi(varargin,'Howfar'))+1};
else
   bnd = ndpr-1; 
end
%
if sum(strcmpi(varargin,'ACtype'))
    ACm           =   varargin{find(strcmpi(varargin,'ACtype'))+1};
    if strcmpi(ACm,'Unbiased')
        ACmfalg=1;
    elseif strcmpi(ACm,'robust')
        ACmfalg=2;
    else
        error('Unknown method!')
    end
else
    ACm='SampleAC';
    ACmfalg=0;
end

%% Check Param

if size(X,2)==1  
    if verbose; disp(['TS1 timeseries transposed!']); end;
    X=X';
end
if size(Y,2)==1
    if verbose; disp(['TS2 timeseries transposed!']); end; 
    Y=Y';
end

assert(length(Y)==length(X),'Length of time series should be equal.')
%% Basic Params
ndpr = length(X);
%j = round(ndpr*Howfar);

%% Autocorrelation Estimation
if  ACmfalg==0 %Biased
    ac_x = autocorr(X,bnd);
    ac_y = autocorr(Y,bnd);
    
elseif ACmfalg==1 %Unbiased
    ac_x = (ndpr./(ndpr-(0:bnd)))'.*autocorr(X,bnd);
    ac_y = (ndpr./(ndpr-(0:bnd)))'.*autocorr(Y,bnd);
    
elseif ACmfalg==2 %Robust -- never mind; too off!
    ac_x(1)=1; ac_y(1)=1;
    for jj=2:bnd+1
        ac_x(jj,1) = madicc(X(1:end-jj),X(jj+1:end));
        ac_y(jj,1) = madicc(Y(1:end-jj),Y(jj+1:end));
    end
end

%% Correction Factor Calculation
wght=(ndpr-(1:bnd));
CF_ind=(ndpr+2.*sum(wght.*(ac_x(2:end).*ac_y(2:end))))./ndpr;

if  CFmethod==1
    CF=CF_ind;
elseif CFmethod==2
    rho = corr(X',Y');
    CF=CF_ind+rho.^2;
    if verbose; disp(['xcorr(X,Y): ' num2str(rho)]); end; 
elseif CFmethod==3        
    wght_cc = (ndpr-abs(-bnd_cc:bnd_cc));
    xycc=crosscorr(X,Y,bnd_cc); %including the \rho_{xy} of course!

    CF_cc=sum(wght_cc.*(xycc.^2))./ndpr;
    CF=CF_ind+CF_cc;
    
    rho=xycc(ceil(length(wght_cc)/2));
    
    if verbose; disp(['-corr(X,Y): ' num2str(rho)]); end; 
end

EDF=ndpr./CF;

if verbose
    disp(['-------------------------------'])
    disp(['-AC Estimation Method: ' ACm]); 
    disp(['-CF Estimation Method: ' CFm]);
    disp(['--How far on ACF?      ' num2str(bnd)]);
    if ~isnan(bnd_cc)
        disp(['--How far on XCorr?    ' num2str(bnd_cc)]);
    end
end

Stat.corrXY     = rho;
Stat.bnd_xcorr  = bnd_cc;
Stat.bnd        = bnd;
Stat.N          = ndpr;
Stat.eN         = EDF;
Stat.BCF        = CF_ind;
Stat.CFm        = CFm;
Stat.ACm        = ACm;


%%%%%%%%%%%%%%%%%%%%%ALTERNATIVE%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Yt1=toeplitz(autocorr(x,length(x)-1)); 
%Yt2=toeplitz(autocorr(y,length(y)-1)); 

%trace(Yt1*Yt2)./length(x) %Effing SLOW
%OR
%sum(sum(Yt1.*Yt2,2))./length(x) %FASTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

