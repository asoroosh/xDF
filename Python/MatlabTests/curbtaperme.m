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