clear

addpath /Users/sorooshafyouni/Home/GitClone/HetBiv/

for i=1:10000
    
    %if ~mod(i,100); disp(i); end;
    
    x1=randn(2000,1); x2=randn(2000,1);
    
    CF(i)=PyperPeterman(x1,x2,'howfar',0.25);
    CF_UB(i)=PyperPeterman(x1,x2,'Autocorr','unbiased','howfar',0.25);
end


figure; 
hold on; box on; 
histogram(CF,70,'Normalization','probability')
ylimh=ylim;
line([mean(CF) mean(CF)],ylimh,'color','r')
histogram(CF_UB,70,'Normalization','probability')
line([mean(CF_UB) mean(CF_UB)],ylimh,'color','b')
legend({'B','meanB','UB','UBmean'})
