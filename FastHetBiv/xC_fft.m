function [xC]=xC_fft(Y,L,varargin)
%[xAC]=xC_fft(Y,L,varargin)
%   Super fast full-lag cross-correlation calculation of multi-dimention 
%   matrices.
%
%   If you need the Pearson's correlation (lag0 xcorr), set lag to 0! Also,
%   diagonal of layer n is the AC of lag n! 
%
%   This function only works for Matlab 2016< and there is no way around
%   it!
%_________________________________________________________________________
% Soroosh Afyouni, NISOx.org, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if size(Y,2)~=L
    Y=Y'; %IxT
end
I = size(Y,1);

if sum(strcmpi(varargin,'lag'))
    mxL = varargin{find(strcmpi(varargin,'lag'))+1};
    if ~mxL; mxL=1; end; %the user wants the Pearson's Correlation!
else
    mxL = L;
end

Y = Y-mean(Y,2);
%only works on >2016 Matlabs, but faster!
% if <2016, use Y=Y-repmat(mean(Y,2),1,L) instead.

nfft = 2^nextpow2(2*L-1);
yfft = fft(Y,nfft,2);

mxLcc = (mxL-1)*2+1;
xC    = zeros(I,I,mxLcc);

[xx,yy]=find(triu(ones(I)));
for i=1:numel(xx)
    xC0         = ifft(yfft(xx(i),:).*conj(yfft(yy(i),:)));
    xC0         = [xC0(end-mxL+2:end),xC0(1:mxL)];
    xC(xx(i),yy(i),:)   = xC0./sqrt(sum(abs(Y(xx(i),:)).^2)*sum(abs(Y(yy(i),:)).^2));
    clear xC0
end

