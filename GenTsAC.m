%function ts=GenTsAC(ar1,ndpr)
% Generates very dumb non-white noises 
% ar1:  The first-order autoregressive coefficient
% ndpr: Number of time points

function ts=GenTsAC(ar1,ndpr)

ts(1)=randn(1); 
for t=2:ndpr 
    ts(t)=ts(t-1)*ar1+randn(1); 
end