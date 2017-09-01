function [BCF,EDF,xAC]=HetBivCalc_fft(Y,L,varargin)
%[BCF,xAC]=HetBivCalc_fft(Y,L)
%
%   Super fast full-lag Bartlett's Correction Factor (BCF) calculation of 
%   multi-dimention matrices. 
%   
%   NB! If input Y is not in IxT form, then it is automatically transposed.
%
%   Dependencies: AC_fft.m & xC_fft.m
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
    if strcmpi(CFm,'BCF')
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
        error('Unknown BCF methods. Choose between BCF, BC0, BCx');
    end
else
   CFmethod=1;
   lagcc=0;
end

xAC      = AC_fft(Y,L);
xAC(:,1) = []; %because we later take care of that little lag0! 
nLg      = L-1;

wgt     = (nLg-(1:nLg));
BCF     = wgt.*xAC*xAC'; %pfff
BCF     = (L+2*(BCF))./L;

if CFmethod==1
    BCF(BCF<1)=1; %ensure that there are no node with BCF larger than 1.
    EDF=L./BCF; 
    return; 
end

if CFmethod==2
    BCF=BCF+corr(Y').^2;
    BCF(BCF<1)=1; %ensure that there are no node with BCF larger than 1.
    EDF=L./BCF; 
end

if CFmethod==3
    xC  = xC_fft(Y,L,'lag',lagcc); %IxIxlagcc
    wgt = (nLg-abs(-(lagcc-1):lagcc-1)); %1xlagcc
    %size(wgt), size(xC)
    xC  = bsxfun(@times,xC.^2,reshape(wgt,1,1,numel(wgt))); %here, wgt is multiplied by each plain of IxI of xC!
    xC  = sum(xC,3);
    xC  = xC+triu(xC,1)';
    BCF = BCF+xC./L;
    
    BCF(BCF<1)=1; %ensure that there are no node with BCF larger than 1.
    
    EDF = L./BCF;
end 


