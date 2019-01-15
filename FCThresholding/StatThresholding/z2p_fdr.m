function [h_fdr,padj_fdr,p_unadj,pval_cv_fdr]=z2p_fdr(z_mat)
%[h_fdr,p_fdr,p,cv_fdr]=z2p_fdr(z_mat)
%
%%%INPUTS:
% where z_mat for naive correction factor is atanh(r).*sqrt(df-3);
%
% Input can be either matrix of NxN or vector of 1xN. 
% If input is 3D, the 3rd dim is considered as number of subjects, and
% the std is adjusted accordingly. 
%
%%%OUTPUTS:
%
%   h_fdr   : flags of whether an edge is sig or not.
%   p_fdr   : adjusted p-values (with fdr, obviousely!).
%   p_unadj : raw pvalues (i.e. unadjusted).
%   pval_cv_fdr : critical pvalues on FDR corrected.
%
%%%DEPENDENCIES:
% fdr_bh.m
%
%%%REFERENCES:
%   Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%   under Serial Correlation
%   Soroosh Afyouni, Stephen M. Smith & Thomas E. Nichols
%   2018
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
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

p_unadj = 2.*normcdf(-abs(z_mat),0,ms_std);

%just to make sure that digonal won't be detected!!
if itsamat; nn=size(p_unadj,1); p_unadj(1:nn+1:end) = 1; end

[Hflag_fdrx,pval_cv_fdr,adj_p] = fdr_bh(p_unadj(idx));

%matrix of adjusted p-values
padj_fdr       = zeros(size(p_unadj));
padj_fdr(idx)  = adj_p;
padj_fdr       = padj_fdr+padj_fdr';

%matrix of adjusted p-values - in binary
h_fdr       = zeros(size(p_unadj));
h_fdr(idx)  = Hflag_fdrx;
h_fdr       = h_fdr+h_fdr';