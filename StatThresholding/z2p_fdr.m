function [h_fdr,p_fdr,p,cv_fdr]=z2p_fdr(z_mat)
%[h_fdr,p_fdr,p,cv_fdr]=z2p_fdr(z_mat)
%
%%%INPUTS:
% where z_mat for naive correction factor is atanh(r).*sqrt(df-3);
%
% Input can be either matrix of NxN or vector of 1xN. 
% If input is 3D, the 3rd dim is considered as number of subjects, and
% the std is adjusted accordingly. 
%
%%%DEPENDENCIES:
% fdr_bh.m
%_________________________________________________________________________
% Soroosh Afyouni, NISOx.org, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

%Notes: Checked it with corrcoeff.m for null case, 1 f thousand times!
itsamat=0;
if size(z_mat,1)==size(z_mat,2) && size(z_mat,1)~=1 %ohhhh! for matrices!
    itsamat=1;
    idx = find(triu(z_mat,1)); %only the upper triangle
elseif size(z_mat,1)==1 || size(z_mat,2)==1
    idx = 1:numel(z_mat); %and here for vectors!
elseif numel(z_mat)==1
    idx=1; % single number
else
    error('Check the input!')
end
%nn = size(z_mat,1);
ns = size(z_mat,3);
if ns>1
    warning('Input is multi-subject, we''re gonna do group-analysis now!')
    ms_std=1./sqrt(ns);
    z_mat=mean(z_mat,3);
else 
    ms_std=1;
end

p = 2*normcdf(-abs(z_mat),0,ms_std);

%just to make sure that digonal won't be detected!!
if itsamat; nn=size(p,1); p(1:nn+1:end) = 1; end

[p_fdrx,cv_fdr,adj_p] = fdr_bh(p(idx));

%matrix of adjusted p-values
p_fdr       = zeros(size(p));
p_fdr(idx)  = adj_p;
p_fdr       = p_fdr+p_fdr';

%matrix of adjusted p-values - in binary
h_fdr       = zeros(size(p));
h_fdr(idx)  = p_fdrx;
h_fdr       = h_fdr+h_fdr';