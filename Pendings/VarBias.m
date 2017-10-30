clear

BCF=@(a,b,N) (N+2.*sum(((N-1):-1:1).*(a(2:end).*b(2:end))));

Nlist=[10 50 100 1000 10000];
%Nlist=5000
rholist=[0:0.2:0.8 0.90];
nRlz=1000;

varscale=((1-rholist.^2).^2)./(1+rholist.^2);

R = zeros(nRlz,numel(Nlist),numel(rholist));
P = zeros(nRlz,numel(Nlist),numel(rholist));
rr=1;
for rho=rholist
    nn=1;
    for N=Nlist
        disp([num2str(rho) '-' num2str(N)])
        for i=1:nRlz 
            ts = mvnrnd([0 0],[1 rho; rho 1],N);
            
            R(i,nn,rr) = corr(ts(:,1),ts(:,2));            
            P(i,nn,rr) = ts(:,1)'*ts(:,2)/N;
            
            %xAC = AC_fft(ts,N);
            %vPt(i,nn,rr) = BCF(xAC(1,:),xAC(2,:),N)/N^2;
            
            [~,v]=HetBivCalc_fft(mvnrnd([0 0],[1 rho; rho 1],N),N,'method','CF0');
            HetBiv(i,nn,rr)=v(1,2);
        end
        nn = nn+1;
    end
    rr = rr+1;
end

disp('%%%% Var(r_xy) vs Var(X`Y)/N %%%%')
stdP=squeeze(var(P)).*repmat(varscale,5,1);
stdR=squeeze(var(R));
disp('%%%%%% BIAS %%%%%%%%%%%%%%%%%%%%%')
(stdR-stdP)./stdR*100


disp('%%%%%% EDOF %%%%%%%%%%%%%%%%%%%%%')
Pedf=1+1./(stdP);
Redf=1+1./(stdR);
%round(squeeze(1+1/var(atanh(R))))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
tP=r.*sqrt((Pedf-2)./(1-rho.^2));



clear R0 P0 R P ts V V0 Rv
 BCF=@(a,b,N) (N+2.*sum(((N-1):-1:1).*(a(2:end).*b(2:end))));
 N=1000; rho=0;
 for i=1:10000
     %ts=mvnrnd([0 0],[1 rho; rho 1],N); 
     ts=corrautocorr([0 0],rho,[0.4 0.2],N)';
     
     P0(i)=(ts(:,1)'*ts(:,2)/N); 
     R0(i)=corr(ts(:,1),ts(:,2));
     
     [~,v]=HetBivCalc_fft(ts,N);
     v=v(1,2);
     
     %xAC=AC_fft(ts,N);
     %v=BCF(xAC(1,:),xAC(2,:),N)./N^2;
     
     V0(i)=v;
     
     vr=(N.*rho.^2./N.^2); % trace(sigma_xy^2)
     vr=vr.*((1-rho.^2).^2)/(1+rho.^2); 
     %V(i)=(v+vr).*((1-rho.^2).^2)/(1+rho.^2);
     V(i)=(v+vr);
     
     clear v vr
 end; 
%Pv=var(P0).*((1-rho^2)^2)/(1+rho^2) 
Rv=var(atanh(R0));
disp(['MC    :' num2str(Rv)])
disp(['Depndt:' num2str(mean(V))])
disp(['Null  :' num2str(mean(V0))])
%(Rv-Pv)/Rv*100
%trace(toeplitz(xAC(1,:))*toeplitz(xAC(2,:)))./N^2


% Pedf=1+1./(Pv);
% tP=rho.*sqrt((Pedf-2)./(1-rho.^2));
% tn=rho.*sqrt((N-2)./(1-rho.^2));
% 
% (tn-tP)/tn*100


%figure; 
%plot(squeeze(var(P))./squeeze(var(atanh(R))))
%figure;
%plot(squeeze(var(P))./squeeze(var(R)))
%figure; 
%hold on; 
%plot(squeeze(mean(vPt)),'r')
%plot(squeeze(var(P)),'b')
