function [CF,EDF,xAC]=HetBivCalc_fft(Y,L,varargin)
%[CF,EDF,xAC]=HetBivCalc_fft(Y,L,varargin)
%
%   Super fast full-lag Bartlett's Correction Factor (BCF) calculation of 
%   multi-dimention matrices. 
%
%%%INPUTS:
%   Y: IxT or TxI BOLD voxel-wise/pacellated time series 
%   L: Length of the time series! Just to check the dimensions are all
%   sorted!
%
%   Optional Inputs:
%       Method: 'BHCF', 'CF0', 'CFx'. [default: BHCF]
%       lag: curbing on cross-covariance structure. [default: 1%]
%
%
%%%OUTPUTS:
%   BCF:  IxI matrix. Correction Factor selected as Method. \hat{N}=N/BCF
%   EDF:  IxI matrix. Effective Degree of Freedom.
%   xAC:  Full-lag Autocorrelation. IxT-1, Becareful about memory! 
%   NB! if the correction factor is below 1, then it is forced to 1.
%%%DEPENDENCIES: 
%   AC_fft.m 
%   xC_fft.m
%
%%%REFERENCES:
%
%   Afyouni & Nichols, 2017
%   Bayley & Hammersley, 1946
%   Richardson & Clifford, 1991
%_________________________________________________________________________
% Soroosh Afyouni, NISOx.org, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if size(Y,2)~=L
    Y=Y'; %IxT
end
I = size(Y,1);

if sum(strcmpi(varargin,'Method'))
   CFm = varargin{find(strcmpi(varargin,'Method'))+1};
    if strcmpi(CFm,'BHCF')
        CFmethod=1;
        lagcc=0;
    elseif strcmpi(CFm,'CF0')
        CFmethod=2;
        lagcc=0;
    elseif strcmpi(CFm,'CFx')
        CFmethod=3;
        if sum(strcmpi(varargin,'lag'))
           lagcc = varargin{find(strcmpi(varargin,'lag'))+1};
        else
           lagcc = round(L*0.01); %because we found that if you keep it short it would work better.
           if ~lagcc; lagcc=1; end;
        end
    else
        error('Unknown CF methods. Choose between BHCF, BC0, BCx');
    end
else
   CFmethod=1;
   lagcc=0;
end

xAC      = AC_fft(Y,L);
xAC(:,1) = []; %because we later take care of that little lag0! 
nLg      = L-1;

wgt = (nLg-(1:nLg));
CF  = wgt.*xAC*xAC'; %pfff
CF  = (L+2*(CF))./L;

if CFmethod==1
    CF(CF<1) = 1; %ensure that there are no node with BCF smaller than 1.
    EDF = L./CF; 
    return; 
end

if CFmethod==2
    CF = CF+corr(Y').^2;
    CF(CF<1) = 1; %ensure that there are no node with BCF smaller than 1.
    EDF = L./CF; 
end

if CFmethod==3
    xC  = xC_fft(Y,L,'lag',lagcc); %IxIxlagcc
    wgt = (nLg-abs(-(lagcc-1):lagcc-1)); %1xlagcc
    %size(wgt), size(xC)
    xC  = bsxfun(@times,xC.^2,reshape(wgt,1,1,numel(wgt))); %here, wgt is multiplied by each plain of IxI of xC!
    xC  = sum(xC,3);
    xC  = xC+triu(xC,1)';
    CF  = CF+xC./L;
    
    CF(CF<1) = 1; %ensure that there are no node with BCF smaller than 1.
    
    EDF = L./CF;
end 


