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
    
    (1+cos([1:M].*pi./M))./2
    
    tt_ts(1:M) = (1+cos([1:M].*pi./M))./2.*acs(1:M);
    %figure; plot(tt_ts); hold on; plot(acs); 
end