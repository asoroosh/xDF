clear

% #########################################################################
% To validate Matlab and Octave versions
% SA, 2021
% #########################################################################

T = 500; 
nn = 2;

% #########################################################################
% Simulate two time series under different conditions ---------------------
% #########################################################################

% Uncomment one of the following sections...

% uncorrelated white noise ################################################
%ts = randn(T,nn)';

% correlated white ########################################################
% If you choose to simulate this scenario, you need xDF added to your path
%
%rhosim = 0.7; 
%ts = corrautocorr([0 0],rhosim,eye(T),T);

% correlated AR1 ##########################################################
% If you choose to simulate this scenario, you need xDF added to your path
%
rhosim = 0.7; 
ARsim  = 0.5; 
ts = corrautocorr([0 0],rhosim,MakeMeCovMat(ARsim,T),T);

ts = ts-mean(ts,2); 


% Just to check Octave and Matlab version are identical 
[VarOctave,StatOctave] = xDF_octave(ts,T);
[VarMatlab,StatMatlab] = xDF(ts,T);

any(any(VarOctave-VarMatlab))
any(any(StatOctave.z-StatMatlab.z))
any(any(StatOctave.p-StatMatlab.p))