clear 

addpath(genpath('/Users/sorooshafyouni/Home/GitClone/HetBiv/'));
save_dir='/Users/sorooshafyouni/Home/GitClone/HetBiv/Pendings/R';

ndp  = 1200;
nRlz = 5000;

lw=1.3;
fs=12;
nbin=20;

postfix='_l';

edf = @(v) (1./v)-1;
c2t = @(C,edf) C.*sqrt((edf-2)./(1-C.^2));

ac_95bnd=[-2/sqrt(ndp) 2/sqrt(ndp)];

MethodList={'Clliford','Dutilleul'};

rho_lis = [  0     0.4   0.8];
Rholev  = {'None','Mid','High'};

%AC_list = {eye(ndp),[0.6:-0.1:0.1]};
AClev   = {  'None',   'Full'};


hf0=figure('position',[50,500,800,600]); 
spc=1;
for rhoIdx=1:3
    for acIdx=[1 2]
        subplot(3,2,spc)
        grid on; hold on; box on; 
        title(['AC: ' AClev{acIdx} '   $\rho$ :' num2str(rho_lis(rhoIdx))],'fontsize',fs,'Interpreter','Latex')
        m_cnt=1;
        for Method=MethodList
            load([save_dir '/MVNR_Sims_' Method{1} '_' AClev{acIdx} 'ac_' Rholev{rhoIdx} 'rho' postfix '.mat'],'V')
            [cnt,num]=HistLine(V,nbin);
            plot(cnt,num,'linewidth',lw)
            
            xlim([0.5 8.5]);
            m_cnt=m_cnt+1;
        end
        spc=spc+1;
        ylabel('Density','fontsize',fs,'Interpreter','Latex')
        xlabel('Correction Factors','fontsize',fs,'Interpreter','Latex')
        legend(MethodList,'fontsize',fs,'Interpreter','Latex')
    end
end
set(hf0,'Color','w')
%export_fig(hf0,['Figs/MatNormSim_CF.pdf'])


hf1=figure('position',[50,500,800,300]); 
hold on; grid on; box on; 
spc=1;
for acIdx=[1 2]
    subplot(1,2,acIdx); 
    box on; hold on; grid on; 
    title(['Autocorrelation: ' AClev{acIdx}])
    load([save_dir '/MVNR_Sims_' MethodList{1} '_' AClev{acIdx} 'ac_' Rholev{2} 'rho' postfix '.mat'],'xAC')
    plot(xAC(1:100:end,:)','color',[.5 .5 .5 .5])
    l0=plot(mean(xAC),'color','r','marker','o');
    l1=line([1 100],[ac_95bnd(1) ac_95bnd(1)],'linewidth',lw,'linestyle','-.');
    line([1 100],[ac_95bnd(2) ac_95bnd(2)],'linewidth',lw,'linestyle','-.');
    %ylim([-.7 .7])
    %xlim([1 100])
    
    legend([l0 l1],{'<ACF>','CI 95%'},'fontsize',fs)
    
    ylabel('Autocorrelation Coefficients','fontsize',fs,'Interpreter','Latex')
    xlabel('Lags','fontsize',fs,'Interpreter','Latex')
end
set(hf1,'Color','w')
%export_fig(hf1,['Figs/MatNormSim_AC.pdf'])

hf2=figure('position',[50,500,800,600]); 
spc=1;
for rhoIdx=1:3
    for acIdx=[1 2]
        subplot(3,2,spc)
        grid on; hold on; box on; 
        title(['AC: ' AClev{acIdx} '   $\rho$ :' num2str(rho_lis(rhoIdx))],'fontsize',fs,'Interpreter','Latex')
        m_cnt=1;
        for Method=MethodList
            load([save_dir '/MVNR_Sims_' Method{1} '_' AClev{acIdx} 'ac_' Rholev{rhoIdx} 'rho' postfix '.mat'],'C')
            [cnt,num]=HistLine(C,nbin);
            plot(cnt,num,'linewidth',lw)
            
            xlim([-1 1]);
            m_cnt=m_cnt+1;
        end
        spc=spc+1;
        ylabel('Density','fontsize',fs,'Interpreter','Latex')
        xlabel('Correction Coefficients','fontsize',fs,'Interpreter','Latex')
        legend(MethodList,'fontsize',fs,'Interpreter','Latex')
    end
end
set(hf2,'Color','w')
%export_fig(hf2,['Figs/MatNormSim_CorrCoeff.pdf'])

% xlimlist=[-6 6;2 22;10 45];
% hf3=figure('position',[50,500,800,600]); 
% spc=1;
% for rhoIdx=1:3
%     for acIdx=[1 2]
%         subplot(3,2,spc)
%         grid on; hold on; box on; 
%         title(['AC: ' AClev{acIdx} '   $\rho$ :' num2str(rho_lis(rhoIdx))],'fontsize',fs,'Interpreter','Latex')
%         m_cnt=1;
%         for Method=MethodList
%             load([save_dir '/MVNR_Sims_' Method{1} '_' AClev{acIdx} 'ac_' Rholev{rhoIdx} 'rho' postfix '.mat'],'C','CF')
%             
%             Cz=atanh(C).*sqrt(ndp./CF);
%             [cnt_z,num_z]=HistLine(Cz,nbin);
%             plot(cnt_z,num_z,'linewidth',lw)
%             %xlim([-1 1]);
%             m_cnt=m_cnt+1;
%         end
%         
%         CzN=atanh(C).*sqrt(ndp-3);
%         [cnt_zN,num_zN]=HistLine(CzN,nbin);
%         plot(cnt_zN,num_zN,'linewidth',lw,'color',[.5 .5 .5])        
%         
%         ylabel('Density','fontsize',fs,'Interpreter','Latex');
%         xlabel('Fisher Z','fontsize',fs,'Interpreter','Latex');
%         
%         if spc==1
%             legend([MethodList,'Naive'],'fontsize',fs,'Interpreter','Latex');
%         end
%         
%         xlim(xlimlist(rhoIdx,:))
%         
%         spc=spc+1;
%     end
% end
% set(hf3,'Color','w')
%export_fig(hf3,['Figs/MatNormSim_FihserZ.pdf'])

xlimlist=[-6 6;2 22;10 45];
hf3=figure('position',[50,500,800,600]); 
spc=1;
for rhoIdx=1:3
    for acIdx=[1 2]
        subplot(3,2,spc)
        grid on; hold on; box on; 
        title(['AC: ' AClev{acIdx} '   $\rho$ :' num2str(rho_lis(rhoIdx))],'fontsize',fs,'Interpreter','Latex')
        m_cnt=1;
        for Method=MethodList
            load([save_dir '/MVNR_Sims_' Method{1} '_' AClev{acIdx} 'ac_' Rholev{rhoIdx} 'rho' postfix '.mat'],'V')
            
            Ct=c2t(C,edf(V));
            [cnt_z,num_z]=HistLine(Ct,nbin);
            plot(cnt_z,num_z,'linewidth',lw)
            %xlim([-1 1]);
            m_cnt=m_cnt+1;
        end
        
        CtN=c2t(C,ndp);
        [cnt_zN,num_zN]=HistLine(CtN,nbin);
        plot(cnt_zN,num_zN,'linewidth',lw,'color',[.5 .5 .5])        
        
        ylabel('Density','fontsize',fs,'Interpreter','Latex');
        xlabel('Student t','fontsize',fs,'Interpreter','Latex');
        
        if spc==1
            legend([MethodList,'Naive'],'fontsize',fs,'Interpreter','Latex');
        end
        
        %xlim(xlimlist(rhoIdx,:))
        
        spc=spc+1;
    end
end
set(hf3,'Color','w')
%export_fig(hf3,['Figs/MatNormSim_FihserZ.pdf'])