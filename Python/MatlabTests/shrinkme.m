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