% function [RtoZcorrection,arone_n,grotZ]=AR1MC_HCPS900_CatTs(TS,nn,r,d,ndpr)
% This function is AR(1) Monte-Carlo estimation of unbiased variance. 
% Copy-pasted from FSLnets toolbox. 
%
% NB! This function was designed for estimation of *global* correction (subjects)
% factors, in contrast to what HetBiv family offers which is a local (edge)
% corrections.
%
% In FSLnets settings:
% TS: 
% nn:   number of nodes
% r:    number of runs
% d:    number of k-space directions
% ndpr: number of time-points
%
%
% RtoZcorrection: Scales your correlation coefficients 
% arone_n:  Median of AR(1) of all nodes within a subject
% grotZ:    Fisher's transformed correlation coefficients of non-white
%           noises
%



function [RtoZcorrection,arone_n,grotZ]=AR1MC_HCPS900_CatTs(TS,nn,r,d,ndpr)

Nruns=r*d;
ndpo=ndpr*Nruns;

arone=[];
for s=1:Nruns
    grot=TS((s-1)*ndpr+1:s*ndpr,:);
    for i=1:nn %SA: calc arone for each node
      g=grot(:,i);  arone=[arone sum(g(1:end-1).*g(2:end))/sum(g.*g)];
    end
    % end
    arone_n(s)=median(arone); %get the median of the nodes AR1
end
arone=median(arone_n); %get the median for all runs

% create null data using the estimated AR(1) coefficient
clear grot*; grotRR=[];
% for s=1:ts.Nsubjects
for i=1:nn 
  grot(1)=randn(1);
  for t=2:ndpo %SA: creat random numbers with ar1 element for rest of the *run*
    grot(t)=grot(t-1)*arone+randn(1);
  end
  grotts(:,i)=grot;
end

grotr=corr(grotts);

grotR=grotr(eye(nn)<1);
grotZ=0.5*log((1+grotR)./(1-grotR));
RtoZcorrection=1/std(grotZ);



