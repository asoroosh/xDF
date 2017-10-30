clear

% tlist=[20 200 600 1200];
% aclist={0,0.4,[0.5 0.4 0.2],[-0.4 0.4]};
% rholist=[0.1 0.5 0.8];
% 
% nRlz=2000;

load('Sim_ARCD_5k_HetBiv2.mat')

fp=figure('position',[50,500,800,600]); hold on; 
for t_cnt=1:numel(tlist)
    

    sp0=subplot(2,2,t_cnt);
    hold on; grid on; box on; 
    
    title(['T= ' num2str(tlist(t_cnt))])
    c=sqrt(1./squeeze(v_c_CloseLags(t_cnt,:,:)));
    mc=sqrt(1./squeeze(v_mc(t_cnt,:,:)));
    
    bias=((mc-c)./mc)*100';
    plot(bias,'marker','o')
    
    ylabel('Bias $\sigma_{xy}, 100\times[\sigma^{MC}-tr(\Sigma_x\Sigma_y)]/\sigma^{MC}]$','Interpreter','latex','fontsize',12)
    xlabel('Autocorrelation Level','Interpreter','latex','fontsize',12)
    
    set(sp0,'Xtick',1:4,'XTickLabel',{'None','AR1','AR4','AR6+negs'},'XTickLabelRotation',45)
    
    ylim([-70 70]); 
    xlim([0.5 4.5])
    line([.5 4.5],[0 0],'color','r')
    
    legend({'\rho=0','\rho=0.1','\rho=0.5','\rho=0.8'},'fontsize',12)
    
end

set(gcf,'Color','w')
%export_fig Sigmaxy_bias.pdf