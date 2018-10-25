function [h_bon,p_bon,p_unadj,alp_bon]=z2p_bon(z_mat)
%[h_bon,p_bon,p,alp_bon]=z2p_bon(z_mat)
%
%%%INPUTS:
% where z_mat for naive correction factor is atanh(r).*sqrt(df-3);
%
% Input can be either matrix of NxN or vector of 1xN. 
% If input is 3D, the 3rd dim is considered as number of subjects, and
% the std is adjusted accordingly. 
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
alp     = .05;
itsamat = 0;
if size(z_mat,1)==size(z_mat,2) && size(z_mat,1)~=1 %ohhhh! for matrices!
    itsamat = 1;
    idx     = find(triu(z_mat,1)); %only the upper triangle
elseif size(z_mat,1)==1 || size(z_mat,2)==1
    idx = 1:numel(z_mat); %and here for vectors!
elseif numel(z_mat)==1
    idx=1; % single number
else
    error('Check the input!')
end
ns = size(z_mat,3);

if ns>1
    warning('Input is multi-subject, we''re gonna do group-analysis now!')
    ms_std = 1./sqrt(ns);
    z_mat  = mean(z_mat,3);
else 
    ms_std = 1;
end

nt      = length(idx); %number of tests
p_unadj       = 2*normcdf(-abs(z_mat),0,ms_std);

%just to make sure that digonal won't be detected!!
if itsamat; nn = size(p_unadj,1); p_unadj(1:nn+1:end) = 1; end

p_bon   = p_unadj.*nt;
alp_bon = alp./nt;

idx_m        = p_bon<alp;
h_bon        = zeros(size(p_unadj));
h_bon(idx_m) = 1;

%pff, are you nuts?! just send back the p-val/fullgraph and threshold on 5%!
% alp_bon=alp./length(idx);
% p_bon=p(idx);
% p_bon(p_bon>alp_bon)=0;
% h_bon=zeros(size(p_bon));
% h_bon(idx)=p_bon;
% h_bon((h_bon.*100e5)>0)=1;
% h_bon=h_bon+h_bon';