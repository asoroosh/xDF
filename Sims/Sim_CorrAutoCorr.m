%
% $Id: dudude$ 
clear


addpath /Users/sorooshafyouni/Home/matlab/Ext/spm12
addpath /Users/sorooshafyouni/Home/GitClone/HetBiv/HetBivs
addpath /Users/sorooshafyouni/Home/GitClone/HetBiv

ndp=500;
rrcor=0.8;

nRlz=1200;
Sigma=[0 rrcor; rrcor 0];
Sigma(eye(size(Sigma))~=0)=1;

for itr=1:nRlz
    if ~mod(itr,100); disp(num2str(itr)); end;
    
    Y=mvnrnd(zeros(1,size(Sigma,1)),Sigma,ndp);
    
    r_w_tmp=corrcoef(Y);
    rw(itr)=r_w_tmp(1,2);
    zh_w_tmp=atanh(r_w_tmp(1,2));

    z_white(itr)=zh_w_tmp;
    
    BCF(itr)=HetBivCalc(Y(:,1)',Y(:,2)','Method','BCF','Howfar',round(0.50*ndp));
    BCx(itr)=HetBivCalc(Y(:,1)',Y(:,2)','Method','BCx','Howfar',round(0.50*ndp),'CCHowfar',round(0.25*ndp));
    BC0(itr)=HetBivCalc(Y(:,1)',Y(:,2)','Method','BC0','Howfar',round(0.50*ndp));
    
    [BCF35_tmp,BCF35A_tmp]=BCF35nArbab_HCPS900_CatTs_TsComb(Y(:,1),Y(:,2),r_w_tmp(1,2),ndp);
    BCF35(itr) = BCF35_tmp;
    BCF35A(itr)= BCF35A_tmp;
    
    z_w_naive(itr)  = zh_w_tmp.*sqrt(ndp-3);
    z_w_bcf(itr)    = zh_w_tmp.*sqrt((ndp./BCF(itr))-3);
    z_w_bc0(itr)    = zh_w_tmp.*sqrt((ndp./BC0(itr))-3);
    z_w_bcx(itr)    = zh_w_tmp.*sqrt((ndp./BCx(itr))-3);
    
    z_w_bcf35(itr)    =  zh_w_tmp.*sqrt((ndp./BCF35(itr))-3);
    z_w_bcf35a(itr)    = zh_w_tmp.*sqrt((ndp./BCF35A(itr))-3);
    
    %BCFw(itr)=HetBiv_ShortBC_Fast(Y(:,1),Y(:,2),ndp,0.15);
    %[BCF35w_tmp,BCF35Aw_tmp]=BCF35nArbab_HCPS900_CatTs_TsComb(Y(:,1),Y(:,2),r_whitetmp(1,2),ndp);
    %BCF35w(itr)=BCF35w_tmp;
    %BCF35Aw(itr)=BCF35Aw_tmp;
    %% Colour them up!
    %ARs=exp(-1*PLc*(1:ndp));

    ARs=[0.8,-0.3,0.05];
    ARs=ARs(round(ARs,1)>0);

    AR1=spm_Q(ARs,ndp,0); 
    ColorK=sqrtm(full(AR1));
    Yc=ColorK*Y;
% 
    %% Validations
    r_c_tmp=corrcoef(Yc);
    rc(itr)=r_c_tmp(1,2);
    zh_c_tmp=atanh(r_c_tmp(1,2));
    
    z_col(itr)=zh_c_tmp;
    
    BCF_c(itr)=HetBivCalc(Yc(:,1)',Yc(:,2)','Method','BCF','Howfar',round(0.50*ndp));
    BCx_c(itr)=HetBivCalc(Yc(:,1)',Yc(:,2)','Method','BCx','Howfar',round(0.50*ndp),'CCHowfar',round(0.25*ndp));
    BC0_c(itr)=HetBivCalc(Yc(:,1)',Yc(:,2)','Method','BC0','Howfar',round(0.50*ndp));    
    
    [BCF35_c_tmp,BCF35A_c_tmp]=BCF35nArbab_HCPS900_CatTs_TsComb(Yc(:,1),Yc(:,2),r_c_tmp(1,2),ndp);
    BCF35_c(itr) = BCF35_c_tmp;
    BCF35A_c(itr)= BCF35A_c_tmp;
    
    
    z_c_naive(itr)  = zh_c_tmp.*sqrt(ndp-3);
    z_c_bcf(itr)    = zh_c_tmp.*sqrt((ndp./BCF_c(itr))-3);
    z_c_bc0(itr)    = zh_c_tmp.*sqrt((ndp./BC0_c(itr))-3);
    z_c_bcx(itr)    = zh_c_tmp.*sqrt((ndp./BCx_c(itr))-3);

    z_c_bcf35(itr)    =  zh_c_tmp.*sqrt((ndp./BCF35_c(itr))-3);
    z_c_bcf35a(itr)    = zh_c_tmp.*sqrt((ndp./BCF35A_c(itr))-3);
%     BCF(itr)=HetBiv_ShortBC_Fast(Yc(:,1),Yc(:,2),ndp,0.15);
%     [BCF35_tmp,BCF35A_tmp]=BCF35nArbab_HCPS900_CatTs_TsComb(Yc(:,1),Yc(:,2),rtmp(1,2),ndp);
%     BCF35(itr)=BCF35_tmp;
%     BCF35A(itr)=BCF35A_tmp;   
    
    clear *tmp*
end

figure; 
hold on
histogram(BC0,50)
histogram(BCx,50)
histogram(BCF,50)
legend({'BC0','BCx','BCF'})

% figure; 
% hold on
% histogram(BCF35w,50)
% histogram(BCF35Aw,50)
% histogram(BCFw,50)
% legend({'BCF35w','BCF35Aw','BCFw'})
% 
% figure; 
% hold on; 
% histogram(r,50)
% histogram(rw,50)
% legend({'r','r-White'})
% 

figure; 
hold on
title('White')
histogram(z_w_naive ,50 ,'normalization','probability')

histogram(z_w_bcf35 ,50 ,'normalization','probability')
histogram(z_w_bcf35a,50 ,'normalization','probability')

histogram(z_w_bcf ,50 ,'normalization','probability')
histogram(z_w_bc0 ,50 ,'normalization','probability')
histogram(z_w_bcx ,50 ,'normalization','probability')
%histogram(z_col.*sqrt(ndp-3),50,'normalization','probability')
legend({'Naive','bcd35','bcf35a','z-BCF','z-BC0','z-BCx'})


figure; 
hold on
title('Coloured')

%histogram(z_w_naive ,50 ,'normalization','probability')

histogram(z_c_bcf35 ,50 ,'normalization','probability')
histogram(z_c_bcf35a,50 ,'normalization','probability')

histogram(z_c_bcf  ,50 ,'normalization','probability')
histogram(z_c_bc0  ,50 ,'normalization','probability')
histogram(z_c_bcx  ,50 ,'normalization','probability')

histogram(z_c_naive,50 ,'normalization','probability')

line(median([z_w_naive;z_w_naive],2),ylim,'color','r')
legend({'bcd35','bcf35a','z-BCF','z-BC0','z-BCx','Naive'})