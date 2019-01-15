clear

V = '/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_124422_OnlyMTS.mat';
load(V)
T = 1200;

[V,Stat] = xDF(mts,T,'truncate','adaptive','verbose');

Stat.z(1:5,1:5)

[h_fdr,padj_fdr,p_unadj,pval_cv_fdr] = z2p_fdr(Stat.z);
numel(find(h_fdr))/2

[h_fdr,padj_fdr,p_unadj,pval_cv_fdr] = z2p_bon(Stat.z);