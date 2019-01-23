clear

V = '/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_124422_OnlyMTS.mat';
load(V)
T = 1200;


[V,Stat] = xDF(mts,T);
[h_bon,padj_bon,p_unadj,pval_cv_bon] = z2p_bon(Stat.z);
numel(find(h_bon))/2

[V,Stat] = xDF(mts,T,'truncate',T/4,'verbose');
[h_bon,padj_bon,p_unadj,pval_cv_bon] = z2p_bon(Stat.z);
numel(find(h_bon))/2

[V,Stat] = xDF(mts,T,'truncate','adaptive','verbose');
[h_bon,padj_bon,p_unadj,pval_cv_bon] = z2p_bon(Stat.z);
numel(find(h_bon))/2

[h_bon,padj_bon,p_unadj,pval_cv_bon] = z2p_fdr(Stat.z);
numel(find(h_bon))/2

[V,Stat] = xDF(mts,T,'taper','tukey',sqrt(T),'verbose');
[h_bon,padj_bon,p_unadj,pval_cv_bon] = z2p_bon(Stat.z);
numel(find(h_bon))/2
