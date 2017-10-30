clear 

addpath(genpath('/Users/sorooshafyouni/Home/GitClone/HetBiv/'));
save_dir='R/';

ndp  = 500;
nRlz = 2000;

MethodList={'Dutilleul','Clliford'};

rho_lis = [  0     0.4   0.8];
Rholev  = {'None','Mid','High'};

Amp=0.6; f=5; atcon=5; phi=0; %N=1200;
t=0:round(1/(ndp/10),5):1;
f=f*2*pi;
y0=Amp*cos(f*t + phi).*exp(-atcon*t);
ACstruct=nonzeros(round(y0,2));

%plot(y)
%[0.6:-0.1:0.1]
%AC_list{1} = ACstruct;
%ACstruct = [0.6,-.4,.2,-.1,0.05];

AC_list = {eye(ndp),ACstruct};
AClev   = {  'None',   'Full'};


for rhoIdx=1:3
    for acIdx=[1 2]
           
            rho = rho_lis(rhoIdx);
            AC  = AC_list{acIdx};
            disp([Method{1} ' - rho:' Rholev{rhoIdx} ' - AC:' AClev{acIdx}])
            ts0=corrautocorr([5 5],rho,AC,ndp,0)';
            c0=corr(ts0); C(i)=c0(1,2);
            
         for Method=MethodList
            xAC=[];
            for i=1:nRlz
                
                xAC0 = AC_fft(ts0,ndp); 
                ac1    = xAC0(1,:); 
                ac2    = xAC0(2,:);
                
                B = eye(ndp);
                if strcmp(Method,'Dutilleul')
                    B   = (eye(ndp)-ones(ndp)/ndp)./ndp;
                end;
                
                v_tmp   = trace(B*toeplitz(ac1)*B*toeplitz(ac2))/(trace(B*toeplitz(ac1))*trace(B*toeplitz(ac2)));
                V(i)=v_tmp;

                xAC=[xAC;xAC0(:,1:100)];

                clear xAC0 *_tmp CF0 ts0 c0 
            end

            save([save_dir '/MVNR_Sims_' Method{1} '_' AClev{acIdx} 'ac_' Rholev{rhoIdx} 'rho_l.mat'],'V','C','xAC','AC','rho','Method','ndp','AC_list','rho_lis')
            clear C V xAC
        end
    end
end