% function CorrLeng=AutoCorrLength(TS,HowFarTs)
% Calculates the Correlation Lengths. i.e. measure how bad the things are
% in terms of autocorrelation. Adapted from:
% 
%       Straatsma, T. P., Berendsen, H. J. C., & Stam, A. J. (2016). 
%       Estimation of statistical errors in molecular simulation calculations,
%       8976(June). http://doi.org/10.1080/00268978600100071
%       
%       Consideration: In the original definition, the correlation length
%       is supposed to consider all lags (see Eq. 4), however, thi function 
%       offer curbing as lots of very far lags, in case of BOLD time
%       series, are merely crappy estimation of zero. 
%
%       HowFarTs should be a value between 0 and 1:
%       HowFarTs=round(1./TS) only takes the first lag
%       HowFarTs=round((TS-1)./TS) consider all lags

function CorrLeng=AutoCorrLength(TS,HowFarTs)

% ndpr=1200;

if size(TS,1)<size(TS,2)
    disp(['timeseries transposed!'])
    TS=TS';
end
ndpr=length(TS);

%Not the best way as the if statement 
%should be coming out of the loops
% size(TS)
% round(ndpr.*HowFarTs)
CorrLeng=sum(autocorr(TS',round(ndpr.*HowFarTs)).^2);
