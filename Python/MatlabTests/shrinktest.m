clear


V = '/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_124422_OnlyMTS.mat';
load(V)
T = 1200;
AC = AC_fft(mts,T);

srnkd_ts = shrinkme(AC(1,:),T);