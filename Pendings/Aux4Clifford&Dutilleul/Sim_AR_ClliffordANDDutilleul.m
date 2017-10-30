clear

addpath(genpath('/Users/sorooshafyouni/Home/GitClone/HetBiv'))

tlist=[20 200 600 1200];
%aclist={0,0.4,[0.8 0.6 0.5 0.3],[0.8 0.6 0.4 -0.3 0.3 -0.2]};
aclist={0,0.4,[0.6 0.4 0.3],[-0.4 0.4]};
rholist=[0 0.1 0.5 0.8];

nRlz=5000;

for t_c=1:numel(tlist)
    for ac_c=1:numel(aclist)
        ac  = aclist{ac_c};
        CovMat0=MakeMeCovMat(ac,tlist(t_c));
        CovMats{t_c,ac_c} = CovMat0;
        for rho_c=1:numel(rholist)
            rho = rholist(rho_c);
            t   = tlist(t_c);
            disp(['t: ' num2str(t) ' ac: ' num2str(ac) 'rho:' num2str(rho)])
                for ii=1:nRlz
                    ts = corrautocorr([0 0],rho,CovMat0,t,1);
                    c_mc00 = corr(ts');
                    c_mc0(ii)  = c_mc00(1,2);
                    clear c_mc00
                    
                    [~,v00]=HetBivCalc_fft(ts,t);
                    v_c0_CloseLags(ii)=v00(1,2);
                    
                    %v_c0(ii)=Clifford_Calc(ts,t);
                    %v_d0(ii)=Dutilleul_Calc(ts,t);
                    
                    %v0_c0(ii)=Clifford_Calc0(ts,t);
                    %v0_d0(ii)=Dutilleul_Calc0(ts,t);
                end

                mcorr(t_c,ac_c,rho_c)= mean(c_mc0);    
                v_mc(t_c,ac_c,rho_c) = var(atanh(c_mc0));
                
                %v_c(t_c,ar_c,rho_c)  = mean(v_c0);
                %v_d(t_c,ar_c,rho_c)  = mean(v_d0); 
                
                v_c_CloseLags(t_c,ac_c,rho_c) = mean(v_c0_CloseLags);
                
                %v0_c(t_c,ar_c,rho_c)  = mean(v0_c0);
                %v0_d(t_c,ar_c,rho_c)  = mean(v0_d0);                 
                
                clear v_c0 v_d0 v0_c0 v0_d0 c_mc0 v_c0_CloseLags
        end
        clear CovMat0
    end
end

save('Sim_ARCD_5k_HetBiv2.mat')