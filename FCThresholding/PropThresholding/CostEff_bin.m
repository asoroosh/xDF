function [CE_val,CE_den,CE_Idx]=CostEff_bin(netmat,densrng,pflag,bflag)
%[CE_val,CE_den,CE_Idx]=CostEff_bin(netmat,densrng,pflag,bflag)
% find cost efficient stat and the corresponding density 
%
% Do you wanna remove the neg values? if yes then pflag = 1;
% Binary procedure? If yes then bflag = 1; 
%
% SA, OX, 2017

if ~exist('bflag','var'); bflag=1; end
if ~exist('pflag','var'); pflag=1; end
if pflag
    netmat(netmat<0)=0; %remove the negs 
end
d_cnt=1;
for d=densrng
    d_net=threshold_proportional(netmat,d);
    if ~bflag
        GE(d_cnt)=efficiency_wei(d_net);
    elseif bflag
        d_net(d_net>0)=1; %binarise
        GE(d_cnt)=efficiency_bin(d_net);
    end
    d_cnt=d_cnt+1;
end
[CE_val,CE_Idx]=max(GE-densrng);
CE_den=densrng(round(mean(CE_Idx)));

