function [V,Z,P,BCF] = CheltonBCF(Y,T)
% The way that Pyper says is right!
% it is fast function i.e. more than 2 time series are accptable for input 
%
% SA, Ox, 2018

if size(Y,1)~=T
    Y=Y'; %IxT
end


ac        = AC_fft(Y,T);
ac(:,[1 end]) = []; %remove the 0lag, we'll take care of it later.

ac = curbtaperme(ac,T-2,round((T-2)/4)); %curb according to Pyper and Peterman

BCF       = 1+2.*(ac*ac');
V         = BCF./T; %is it?
Z         = atanh(corr(Y)).*sqrt(T./BCF-3);
P         = 2 * normcdf(-abs(Z));

end

function ct_ts=curbtaperme(ts,T,M)
% Curb the autocorrelations, according to Anderson 1984
% multi-dimensional, and therefore is fine!
%SA, Ox, 2018
    if ~sum(ismember(size(ts),T)); error('There is something wrong, mate!'); end
    if size(ts,2) ~= T; ts = ts'; end
    
    M          = round(M);
    msk        = zeros(size(ts));
    msk(:,1:M) = 1;
    ct_ts      = msk.*ts;
end