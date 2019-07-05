function [VarHatRho,Stat]=xDF_seedFC(ts_img,ts_sd,T,varargin)
% [VarRho,Stat]=xDF(ts,T,varargin)
% Estimates variance of Pearson's correlations for non-white time series;
%   - Exploites matrix operations & fft for quick estimation of multiple
%   time series. 
%   - Returns variance, adjusted z-scores and uncorrected p-values
%
%
%%%INPUTS:
%   ts_img : Time series of an image as a 2D matrix. 
%   ts_sd : Time series of a seed as a 1D matrix.
%   T : Number of data-points. Just for the sake of sanity checks
%
%   Optionals:
%   'taper'    : uses a tapering method to denoise AC functions
%                Tapering options are:
%                'taper','tukey',M : Single Tukey tapering with cut-off M (as in Woolrich et al 2001)
%
%   'truncate' :
%                'truncate',M: Uses truncations until arbitrary lag M (as in Anderson 1983)
%                'truncate','adaptive', : Uses confidence interval of ACF &
%                remove zero AC coefficients (as suggested in Afyouni, Smith & Nichols 2018)
%
%   'TVOff' : if an estimate exceeed the theoritical variance of a white
%   noise then it curbs the estimate back to (1-rho^2)^2/T. If you want it
%   off, trigger 'TVOff' [default : 'TVOn']
%
%   'verbose': reports some warnings and logs [default: Off]
%
%%%OUTPUTS:
%   VarRho : Variance of rho_ts a 2D matrix of size IxI with diagonal set to zero
%   Stat : is a structure, comprised of:
%        Stat.p: IxI p-values (un-adjusted)
%        Stat.z: IxI z-scores (Fisher transformed & adjusted for AC)
%
%        You probably won't care about the below info:
%        Stat.W2S: IxI of where shrinking has curbed the ACF
%        Stat.TV: Theoritical variance under x & y are i.i.d; (1-rho^2)^2 
%        Stat.EVE: Index of (i,j) edges of which their variance exceeded
%        the theoritical var. 
%
%%%DEPENDECIES:
%   AC_fft.m : estimates ACF super quick via FFT
%   xC_fft.m : estimates cross corr functions super quick via FFT
%
%%%REFERENCES:
% Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. 
% "Effective degrees of freedom of the Pearson's correlation coefficient 
% under autocorrelation." NeuroImage (2019).
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2019
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

    if  size(ts_img,2) ~= T || size(ts_sd,2) ~= T %makes sure dimensions are sound
        warning('Make sure time series are in IxT and are equally sized')
    end
    
    W2S = []; Stat.EVE=[]; TVflag = 1; verbose = 0;
    
    %ts0 : is the image
    %ts1 : is the target time series
    
    nv  = size(ts_img,1); % number of voxels
    
    ts_sd  = ts_sd./std(ts_sd,[],2);
    ts_img = ts_img./std(ts_img,[],2);
    %Corr----------------------------------------------------------------------
    rho   = corr(ts_sd',ts_img')'; %rho is a Ix1 matrix 
    %Autocorr------------------------------------------------------------------
    % AC of the image
    ac_img = AC_fft(ts_img,T); 
    ac_img = ac_img(:,2:T-1); %The last element of ACF is rubbish, the first one is 1, so why bother?!
    % AC of the seed
    ac_sd = AC_fft(ts_sd,T); 
    ac_sd = ac_sd(:,2:T-1); %The last element of ACF is rubbish, the first one is 1, so why bother?!    
    
    nLg  = T-2;         
    %Cross-corr---------------------------------------------------------------- 
    xcf_img  = xC_fft_sd(ts_img,ts_sd,T);
    xc_n = flip(xcf_img(:,2:T-1),2); %negative-lag xcorrs      
    xc_p = xcf_img(:,T+1:end-1);     %positive-lag xcorrs    
    
    %----MEMORY SAVE----
    clear ts_* 
    %-------------------
    
    if sum(strcmpi(varargin,'verbose')); verbose = 1; end     
    
    if sum(strcmpi(varargin,'TVOff'))
        if verbose; disp('Variance Curbing is OFF'); end;
        TVflag = 0; 
    else
        if verbose; disp('Variance Curbing is ON'); end;
    end    
    
    if sum(strcmpi(varargin,'taper'))
        mth = varargin{find(strcmpi(varargin,'taper'))+1};
        if strcmpi(mth,'tukey')
    %Tukey Tappering----------------------------------------
            if nargin<5
                error('you MUST set a tukey factor.');
            elseif ~isnumeric(varargin{find(strcmpi(varargin,'tukey'))+1})
                error('you MUST set a tukey factor.');
            end
            M = round(varargin{find(strcmpi(varargin,'tukey'))+1}); %reads the tukey tapering upper lim (i.e. M; Woolrich et al 2001)
            if verbose; disp(['--Tapering with Tukey of M=' num2str(M)]); end; 
            ac_sd = tukeytaperme(ac_sd,nLg,M);
            for in=1:nv    
                ac_img(in,:) = tukeytaperme(ac_img(in,:),nLg,M);
                xc_n(in,:)   = tukeytaperme(squeeze(xc_n(in,:)),nLg,M);
                xc_p(in,:)   = tukeytaperme(squeeze(xc_p(in,:)),nLg,M);
            end
        else
            error('You can only choose Tukey as tapering option.')
        end
        
    %Truncation------------------------------------------------
     elseif sum(strcmpi(varargin,'truncate'))
         mth = varargin{find(strcmpi(varargin,'truncate'))+1};
        if strcmpi(mth,'adaptive')
            if verbose; disp(['--Adaptive truncation.']); end; 
            for in=1:nv
                    W2S(in) = max([FindBreakPoint(ac_img(in,:),nLg) FindBreakPoint(ac_sd,nLg)]);
            end
            
            ac_sd= shrinkme(ac_sd,nLg);
            for in=1:nv    
                ac_img(in,:)= shrinkme(ac_img(in,:),nLg);
                xc_n(in,:) = curbtaperme(squeeze(xc_n(in,:)),nLg,W2S(in));
                xc_p(in,:) = curbtaperme(squeeze(xc_p(in,:)),nLg,W2S(in));
            end  
    %Curbing------------------------------------------------
        elseif isnumeric(mth)
            M = mth; %the curbing factor
            if verbose; disp(['--Truncation with M=' num2str(M)]); end; 
            ac_sd = curbtaperme(ac_sd,nLg,M);
            for in=1:nv    
                ac_img(in,:) = curbtaperme(ac_img(in,:),nLg,M);
                xc_n(in,:)   = curbtaperme(squeeze(xc_n(in,:))',nLg,M);
                xc_p(in,:)   = curbtaperme(squeeze(xc_p(in,:))',nLg,M);
            end
    %--------------------------------------------------------------------------
        else
            error('Available options are adapative truncation or arbitrary truncation.')
        end
    end
%Crazy, eh?
wgt   = (nLg:-1:1);
wgtm2 = repmat(wgt,[nv,1]);
Tp      = T-1;

 VarHatRho = (Tp*(1-rho.^2).^2 ...
     +   rho.^2 .* sum(wgtm2 .* (SumMat(ac_img.^2,ac_sd.^2,nLg)  +  xc_p.^2 + xc_n.^2),2)...      % 1 2 4
     -   2.*rho .* sum(wgtm2 .* (SumMat(ac_img,ac_sd,nLg)    .* (xc_p    + xc_n))  ,2)...         % 5 6 7 8
     +   2      .* sum(wgtm2 .* (ProdMat(ac_img,ac_sd,nLg)    + (xc_p   .* xc_n))  ,2))./(T^2);   % 3 9 
 
%----MEMORY SAVE----
clear wgtm2 xc_* ac_* 
%-------------------

%Keep your wit about you!
TV = (1-rho.^2).^2./T;
if sum(sum(VarHatRho < TV)) && TVflag
    % Considering that the variance can *only* get larger in presence of autocorrelation.  
    idx_ex       = find(VarHatRho < TV);
    VarHatRho(idx_ex) = TV(idx_ex);
    if verbose; disp([num2str((numel(idx_ex)-nv)/2) ' edges had variance smaller than the textbook variance!']); end;
    Stat.EVE = idx_ex;
end  
Stat.TV    = TV;



%Our turf--------------------------------
rf      = atanh(rho);
sf      = VarHatRho./((1-rho.^2).^2);    %delta method; make sure the N is correct! So they cancel out.
rzf     = rf./sqrt(sf);
f_pval  = 2 .* normcdf(-abs(rzf));  %both tails

%Stat.stable_z = sf;
Stat.z = rzf;
Stat.p = f_pval;

%-------Stat-----------------------------
Stat.W2S = W2S;

end
%--------------------------------------------------------------------------          

%--------------------------------------------------------------------------
function [SM0] = SumMat(Y_img,Y_sd,T)
    %hopefully faster than Matlab's dumb sum of sum!
    %SA, Ox, 2018

    %becareful with this function, this is really tricky to use!
    if ~sum(ismember(size(Y_img),T)) || ~sum(ismember(size(Y_sd),T))
        error('There is something wrong, mate!'); 
    end
    
%     if size(Y_img,1) ~= T || size(Y_sd,1) ~= T
%         error('There is problem with the dimentions.') 
%     end
    
    Y_img = Y_img';
    Y_sd  = Y_sd';
        
    nv  = size(Y_img,2);
    SM0 = zeros(nv,T);

    for i=1:nv        
        SM0(i,:) = (Y_sd+Y_img(:,i));
    end
end

function [SM0] = ProdMat(Y_img,Y_sd,T)
    %hopefully faster than Matlab's dumb sum of sum!
    %SA, Ox, 2018

    %becareful with this function, this is really tricky to use!
    if ~sum(ismember(size(Y_img),T)) || ~sum(ismember(size(Y_sd),T))
        error('There is something wrong, mate!'); 
    end
    
%     if size(Y_img,1) ~= T || size(Y_sd,1) ~= T
%         error('There is problem with the dimentions.') 
%     end
    
    Y_img = Y_img';
    Y_sd  = Y_sd';    
        
    nv  = size(Y_img,2);
    SM0 = zeros(nv,T);
    for i=1:nv
        SM0(i,:) = (Y_sd.*Y_img(:,i));
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
function ct_ts = curbtaperme(acs,T,M)
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
function tt_ts = tukeytaperme(acs,T,M)
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
%--------------------------------------------------------------------------
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
%
%   Stand alone version is in .../xDF/Aux/
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


%--------------------------------------------------------------------------
function [xC,lidx]=xC_fft_sd(ts_img,ts_sd,T)
% Seed based xCF function 
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2017
% srafyouni@gmail.com
% fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if size(ts_sd,2)~=T || size(ts_img,2)~=T 
    error('Make sure the inputs are all IxT matrices/vectors')
end
nv = size(ts_img,1); % size of images

mxL = T;

ts_img = ts_img-mean(ts_img,2);
ts_sd  = ts_sd-mean(ts_sd);
%only works on >2016 Matlabs, but faster!
% if <2016, use Y=Y-repmat(mean(Y,2),1,L) instead.

nfft  = 2^nextpow2(2*T-1); %zero-pad a hell out
fft_img  = fft(ts_img,nfft,2);
fft_sd  = fft(ts_sd,nfft,2);

mxLcc = (mxL-1)*2+1;
xC    = zeros(nv,mxLcc);

for i = 1:nv
    xC0       = ifft(fft_sd.*conj(fft_img(i,:)));
    xC0       = fliplr([xC0(end-mxL+2:end),xC0(1:mxL)]);    
    xC(i,:)   = xC0./sqrt(sum(abs(ts_sd).^2)*sum(abs(ts_img(i,:)).^2));
    clear xC0
end

lidx = [-(mxL-1) : (mxL-1)];

end


