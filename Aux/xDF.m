function [ASAt,Stat]=xDF(ts,T,varargin)
% Estimates the Monster Equation!
%
%%%INPUTS:
%   ts: Time series as a 2D matrix. 
%   T : Number of data-points. Just for the sake of sanity checks
%
%   Optionals:
%   'taper' : uses a tapering method to denoise AC functions
%   'TVOff' : if an estimate exceeed the theoritical variance of a white
%   noise then it curbs the estimate back to (1-rho^2)^2/T. If you want it
%   off, trigger 'TVOff'
%%%OUTPUTS:
%   ASAt : Variance of rho_ts a 2D matrix of size IxI with diagonal set to zero
%   Stat : is a structure, comprised of:
%        Stat.p.f_Pval: IxI p-values (un-adjusted)
%        Stat.z.rzf: IxI z-scores (Fisher transformed & adjusted for AC)
%
%        You probably won't care about the below info:
%        Stat.W2S: IxI of where shrinking has curbed the ACF
%        Stat.TV: Theoritical variance under x & y are i.i.d; (1-rho^2)^2 
%        Stat.EVE: Index of (i,j) edges of which their variance exceeded
%        the theoritical var. 
%%%DEPENDECIES:
%   AC_fft.m : estimates ACF super quick via FFT
%   xC_fft.m : estimates cross corr functions super quick via FFT
%
%%%REFERENCES:
%   Variance of Pearson's correlations under serial-correlations
%   Soroosh Afyouni & Thomas E. Nichols
%   2018
%   University of Oxford
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

    if  size(ts,2) ~= T %makes sure dimensions are sound
        ts = ts';
        warning('Oi!')
    end
    
    W2S = []; Stat.EVE=[]; TVflag = 1; verbose = 0;
    
    nn  = size(ts,1);
    ts  = ts./std(ts,[],2); %standardise
    %Corr----------------------------------------------------------------------
    rho   = corr(ts');
    rho(1:nn+1:end) = 0;
    %Autocorr------------------------------------------------------------------
    [ac] = AC_fft(ts,T); 
    ac   = ac(:,2:T-1); %The last element of ACF is rubbish, the first one is 1, so why bother?!
    nLg  = T-2;         
    %Cross-corr---------------------------------------------------------------- 
    xcf = xC_fft(ts,T);
    xc_n      = flip(xcf(:,:,2:T-1),3); %positive-lag xcorrs
    xc_p      = xcf(:,:,T+1:end-1); %negative-lag xcorrs
    
    %----MEMORY SAVE----
    clear ts 
    %-------------------
    
    if sum(strcmpi(varargin,'TVOff'));   TVflag = 0; end    
    if sum(strcmpi(varargin,'verbose')); verbose = 1; end 
    
    if sum(strcmpi(varargin,'taper'))
        mth = varargin{find(strcmpi(varargin,'taper'))+1};
        if strcmpi(mth,'tukey')
    %Tukey Tappering----------------------------------------
            M = round(varargin{find(strcmpi(varargin,'tukey'))+1}); %reads the tukey tapering upper lim (i.e. M; Woolrich et al 2001)
            if isempty(M); error('you MUST set a tukey factor.'); end;
            for in=1:nn    
                ac(in,:) = tukeytaperme(ac(in,:),nLg,M);
                for jn=1:nn 
                    xc_n(in,jn,:) = tukeytaperme(squeeze(xc_n(in,jn,:)),nLg,M);
                    xc_p(in,jn,:) = tukeytaperme(squeeze(xc_p(in,jn,:)),nLg,M);
                end
            end
            
    %Shrinking------------------------------------------------
        elseif strcmpi(mth,'shrink')
             for in=1:nn
                for jn=1:nn
                    W2S(in,jn) = max([FindBreakPoint(ac(in,:),nLg) FindBreakPoint(ac(jn,:),nLg)]);
                end
            end

            for in=1:nn    
                ac(in,:)= shrinkme(ac(in,:),nLg);
                for jn=1:nn 
                    xc_n(in,jn,:) = curbtaperme(squeeze(xc_n(in,jn,:)),nLg,W2S(in,jn));
                    xc_p(in,jn,:) = curbtaperme(squeeze(xc_p(in,jn,:)),nLg,W2S(in,jn));
                end

            end  
    %Curbing------------------------------------------------
        elseif strcmpi(mth,'curb')
            M = round(varargin{find(strcmpi(varargin,'curb'))+1}); %the curbing factor
            for in=1:nn    
                ac(in,:) = curbtaperme(ac(in,:),nLg,M);
                for jn=1:nn 
                    xc_n(in,jn,:) = curbtaperme(squeeze(xc_n(in,jn,:))',nLg,M);
                    xc_p(in,jn,:) = curbtaperme(squeeze(xc_p(in,jn,:))',nLg,M);
                end
            end
    %--------------------------------------------------------------------------
        else
            error('choose one of these; shrink | tukey | curb as tapering option.')
        end
    end
   
    
    
%Crazy, eh?
wgt     = (nLg:-1:1);
wgtm3   = reshape(repmat((repmat(wgt,[nn,1])),[nn,1]),[nn,nn,numel(wgt)]); %this is shit, eats all the memory!
Tp      = T-1;

 ASAt = (Tp*(1-rho.^2).^2 ...
     +   rho.^2 .* sum(wgtm3 .* (SumMat(ac.^2,nLg)  +  xc_p.^2 + xc_n.^2),3)...         %1 2 4
     -   2.*rho .* sum(wgtm3 .* (SumMat(ac,nLg)    .* (xc_p    + xc_n))  ,3)...         % 5 6 7 8
     +   2      .* sum(wgtm3 .* (ProdMat(ac,nLg)    + (xc_p   .* xc_n))  ,3))./(T^2);   % 3 9 
 
%this part if from PearCorrVarEst.m; when you assume the xCORR are symm! 
% wgtm2   = repmat(wgt,[nn,1]);
% ASAt = [Tp                  .* (1-rho.^2).^2 ...
%        + rho.^2             .* sum(SumMat((wgtm2.*ac.^2),nLg),3) ...                % 1    -- AC 
%        + 2 .* wgt           .* ac*ac'...                                            % 5    -- AC
%        + 2 .* (rho.^2 + 1)  .* sum((wgtm3.*xc_n.*xc_p),3) ...                       %2 & 3 -- XC
%        - 2 .* rho           .* sum(wgtm3.*SumMat(ac,nLg).*(xc_n+xc_p),3)]./(T.^2);  %4 -- This this the only term which we can't seperate the AC and XC! 
%----MEMORY SAVE----
clear wgtm3 xc_* ac 
%-------------------

%Keep your wit about you!
TV = (1-rho.^2).^2./T;
if sum(sum(ASAt < TV)) && TVflag
    % Considering that the variance can *only* get larger in presence of autocorrelation.  
    idx_ex       = find(ASAt < TV);
    ASAt(idx_ex) = TV(idx_ex);
    if verbose; disp([num2str(numel(idx_ex)-nn) ' edges had variance smaller than the textbook variance!']); end;
    [x_tmp,y_tmp]=ind2sub([nn nn],idx_ex);
    Stat.EVE = [x_tmp,y_tmp];
end  
Stat.TV    = TV;

% diagonal is rubbish;
ASAt(1:nn+1:end) = 0;

%------- Test Stat-----------------------
%Pearson's turf -- We don't really wanna go there, eh?
%rz      = rho./sqrt((ASAt));     %abs(ASAt), because it is possible to get negative ASAt!
%r_pval  = 2 * normcdf(-abs(rz)); %both tails
%r_pval(1:nn+1:end) = 0;          %NaN screws up everything, so get rid of the diag, but becareful here. 

%Our turf--------------------------------
rf      = atanh(rho);
sf      = ASAt./((1-rho.^2).^2);    %delta method; make sure the N is correct! So they cancel out.
rzf     = rf./sqrt(sf);
rzf(1:nn+1:end) = 0;
f_pval  = 2 .* normcdf(-abs(rzf));  %both tails
f_pval(1:nn+1:end) = 0;             %NaN screws up everything, so get rid of the diag, but becareful here. 

Stat.z.rzf = rzf;
Stat.p.f_Pval = f_pval;

%Fisher's turf---------------------------
% rf          = atanh(rho);
% edf         = 1./ASAt;                          %Effective Degrees of Freedom
% rzfish      = rf.*sqrt(edf-3);
% rzfish(1:nn+1:end)  = 0;
% f_pval_fish         = 2 .* normcdf(-abs(rzfish));  %both tails
% f_pval_fish(1:nn+1:end) = 0;                %NaN screws up everything, so get rid of the diag, but becareful here. 
% 
% Stat.z.rzfish = rzfish;
% Stat.p.f_PvalFish = f_pval_fish;
%-------Stat-----------------------------
Stat.W2S = W2S;

end
%--------------------------------------------------------------------------          

%--------------------------------------------------------------------------
function [SM0] = SumMat(Y0,T)
    %hopefully faster than Matlab's dumb sum of sum!
    %SA, Ox, 2018

    %becareful with this function, this is really tricky to use!
    if ~sum(ismember(size(Y0),T)); error('There is something wrong, mate!'); end
    if size(Y0,1) ~= T; Y0 = Y0'; end

    nn  = size(Y0,2);
    idx = find(triu(ones(nn),1))';
    SM0 = zeros(nn,nn,T);
    for i=idx
        [x,y]      = ind2sub(nn,i);
        SM0(x,y,:) = (Y0(:,x)+Y0(:,y));
        SM0(y,x,:) = (Y0(:,y)+Y0(:,x));
    end
end

function [SM0] = ProdMat(Y0,T)
    %hopefully faster than Matlab's dumb sum of sum!
    %SA, Ox, 2018

    %becareful with this function, this is really tricky to use!
    if ~sum(ismember(size(Y0),T)); error('There is something wrong, mate!'); end
    if size(Y0,1) ~= T; Y0 = Y0'; end

    nn  = size(Y0,2);
    idx = find(triu(ones(nn),1))';
    SM0 = zeros(nn,nn,T);
    for i=idx
        [x,y]      = ind2sub(nn,i);
        SM0(x,y,:) = (Y0(:,x).*Y0(:,y));
        SM0(y,x,:) = (Y0(:,y).*Y0(:,x));
    end
end

%--------------------------------------------------------------------------

function srnkd_ts=shrinkme(acs,T)
%Shrinks the *early* bucnhes of autocorr coefficients beyond the CI.
%Yo! this should be transformed to the matrix form, those fors at the top
%are bleak!
%
%SA, Ox, 2018
    if ~sum(ismember(size(acs),T)); error('There is something wrong, mate!'); end
    if size(acs,2) ~= T; acs = acs'; end
    %bnd = (sqrt(2)*erfinv(0.95))./sqrt(T);
    %idx = find(abs(acs)>bnd);
    %isit       = abs(acs)>bnd & (1:T);
    %where2stop = find(isit==0); %finds the break point -- intercept 
    %BE CAREFUL:
    %this here is different from Toeplitz ME version, because here we don't
    %have the 0lag, but in that setting the 0lag is there and it is the
    %diag. 
    where2stop = FindBreakPoint(acs,T);
    
    if ~where2stop %if there was nothing above the CI...
        srnkd_ts = zeros(1,T);
    else
        srnkd_ts   = curbtaperme(acs,T,where2stop);
    end
end
%--------------------------------------------------------------------------
function where2stop = FindBreakPoint(acs,T)
% this finds the breaking points for shrinking the AC. 
% Nothing serious, just might help with speed...
% SA, Ox, 2018
    if ~sum(ismember(size(acs),T)); error('There is something wrong, mate!'); end
    if size(acs,2) ~= T; acs = acs'; end
    
    bnd        = (sqrt(2)*erfinv(0.95))./sqrt(T);
    %idx        = find(abs(acs)>bnd);
    isit       = abs(acs)>bnd & (1:T);
    where2stop = find(isit==0); %finds the break point -- intercept 
    
    if where2stop(1)==1
        where2stop = 0; 
    else
        where2stop = where2stop(1)-1; 
    end;
end
%--------------------------------------------------------------------------
function ct_ts=curbtaperme(acs,T,M)
% Curb the autocorrelations, according to Anderson 1984
% multi-dimensional, and therefore is fine!
%SA, Ox, 2018
    if ~sum(ismember(size(acs),T)); error('There is something wrong, mate!'); end
    if size(acs,2) ~= T; acs = acs'; end
    
    M          = round(M);
    msk        = zeros(size(acs));
    msk(:,1:M) = 1;
    ct_ts      = msk.*acs;
end
%--------------------------------------------------------------------------
function tt_ts=tukeytaperme(acs,T,M)
%performs Single Tukey Tapering for given length of window, M, and initial
%value, intv. intv should only be used on crosscorrelation matrices.
%
%NB! There used to be initialisation parameters here before, intv. I 
%remoeved it because we now start with the second elements of the ACF anyways. 
%
%SA, Ox, 2018
    if ~sum(ismember(size(acs),T)); error('There is something wrong, mate!'); end
    if size(acs,2) ~= T; acs = acs'; end
    %if ~exist('intv','var'); intv = 1; warning('Oi!'); end;
    M          = round(M);
    tt_ts      = zeros(size(acs));
    %tt_ts(:,1) = intv;
    tt_ts(1:M) = (1+cos([1:M].*pi./M))./2.*acs(1:M);
    %figure; plot(tt_ts); hold on; plot(acs); 
end

function [xAC,CI,ACOV]=AC_fft(Y,L,varargin)
%[xAC]=AC_fft(Y,T,varargin)
% Super fast full-lag AC calculation of multi-dimention matrices. The
% function exploits fft to estimate the autocorrelations. 
%
%%%%INPUTS
%   Y:      A matrix of size IxT comprised of I time series of T length.
%   L:      Time series length
%   
%   To get a double sided AC, add 'two-sided' as an input.
%
%   NB! The function includes the 0lag zero (i.e. 1) to the output. 
%%%%OUTPUTS
%   xAC:    IxT-1 matrix of full-lag autocorrelations. If 'two-sided'
%           switched, then the matrix is Ix2(T-1).
%   CI :    95% Confidence Intervals of AFC.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if size(Y,2)~=L
    Y=Y';
end

if size(Y,2)~=L
    error('Use IxT or TxI input form.')
end

Y=Y-mean(Y,2); 
%only works on >2016 Matlabs, but faster!
% if <2016, use Y=Y-repmat(mean(Y,2),1,L) instead.

nfft    = 2.^nextpow2(2*L-1); %zero-pad the hell out!
yfft    = fft(Y,nfft,2); %be careful with the dimensions

ACOV = real(ifft(yfft.*conj(yfft),[],2));

ACOV = ACOV(:,1:L);

xAC  = ACOV./sum(abs(Y).^2,2); %normalise the COVs

if sum(strcmpi(varargin,'two-sided')) %two sided is just mirrored, AC func is symmetric
   xAC  = [xAC(:,end-L+2:end) ; xAC(:,1:L)];
else
    xAC  = xAC(:,1:L);
end

bnd=(sqrt(2)*erfinv(0.95))./sqrt(L); %assumes normality for AC
CI=[-bnd bnd];

end

function [xC,lidx]=xC_fft(Y,T,varargin)
%[xAC]=xC_fft(Y,T,varargin)
%   Super fast full-lag cross-correlation calculation of multi-dimensional 
%   matrices.
%
%%%%%  INPUTS:
%   Y:      A matrix of size IxT comprised of I time series of T length.
%   L:      Time series length
%
%%%%%  OUTPUTS:
%   xC   :  Is a 3D matrix of IxIxT: 
%
%           1) On the diagolans, you get the *auto*correlation. Although note
%              that the result is a symmetric ACF (negatives and positives)
%              autocorrelation lags.
%
%           2) Off diagonals are *cross*correlations between a pair
%           3) The identity IxIx1 is a correlation matrix (i.e. lag-0 structures).
%
%   lidx :  I a vector of lag indexes. Each 1x1xT structure follows these
%           lag indexes.
%%%%%% NOTES:
%   If you need the Pearson's correlation (lag0 xcorr), set lag to 0! Also,
%   diagonal of layer n is the AC of lag n! 
%
%   This function only works for Matlab 2016< and there is no way around
%   it!
%
%   For only a pair time series, this is slower than crosscorr. Only use 
%   this function if you have a lager number of time series 
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if size(Y,2)~=T
    Y=Y'; %IxT
end
I = size(Y,1);

if sum(strcmpi(varargin,'lag'))
    mxL = varargin{find(strcmpi(varargin,'lag'))+1};
    if ~mxL; mxL=1; end; %the user wants the Pearson's Correlation!
else
    mxL = T;
end

Y = Y-mean(Y,2);
%only works on >2016 Matlabs, but faster!
% if <2016, use Y=Y-repmat(mean(Y,2),1,L) instead.

nfft  = 2^nextpow2(2*T-1); %zero-pad me
yfft  = fft(Y,nfft,2);

mxLcc = (mxL-1)*2+1;
xC    = zeros(I,I,mxLcc);

[xx,yy]=find(triu(ones(I),1));
for i=1:numel(xx)
    xC0         = ifft(yfft(xx(i),:).*conj(yfft(yy(i),:)));
    xC0         = fliplr([xC0(end-mxL+2:end),xC0(1:mxL)]);
    xC(xx(i),yy(i),:)   = xC0./sqrt(sum(abs(Y(xx(i),:)).^2)*sum(abs(Y(yy(i),:)).^2));
    clear xC0
end

xC = xC + permute(xC,[2 1 3]);

lidx = [-(mxL-1) : (mxL-1)];

end

