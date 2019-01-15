function where2stop = FindBreakPoint(acs,T)
% this finds the breaking points for shrinking the AC. 
% Nothing serious, just might help with speed...
% SA, Ox, 2018
    if ~sum(ismember(size(acs),T)); error('There is something wrong, mate!'); end
    if size(acs,2) ~= T; acs = acs'; end
    
    bnd        = (sqrt(2)*erfinv(0.95))./sqrt(T)
    %idx        = find(abs(acs)>bnd);
    isit       = abs(acs)>bnd & (1:T);
    where2stop = find(isit==0); %finds the break point -- intercept 
    
    if where2stop(1)==1
        where2stop = 0; 
    else
        where2stop = where2stop(1)-1; 
    end;
end