%% Introduction 
% This file is the main script to run full PHY simulation
% shown in Fig.3 and Fig.6 in the IEEE TCOM paper:
% "Efficient PHY Layer Abstraction for Fast Simulations in Complex 
% System Environments"
% The full PHY simulation configurations in this script cover: 
% 1) 20MHz OFDM/OFDMA MIMO/MU-MIMO system
% 2) No interference 
%% Full PHY simulation setup
clear all
tStart = tic;
mcs = [4]; % Vector of MCS to simulate between 0 and 11
numTxRx = [4 2]; % Matrix of MIMO schemes, each row is [numTx numRx]
chan = "Model-D"; % String array of delay profiles to simulate
userIdx = 1; % User of investigation
ruIdx = 1; % RU of investigation
Nsts = 2; % Number of space-time streams
maxnumberrors = 40*1e3;  % The maximum number of packet errors at an SNR point
maxNumPackets = 40*1e3; % The maximum number of packets at an SNR point
% maxnumberrors = 5*1e1;  % The maximum number of packet errors at an SNR point
% maxNumPackets = 5*1e2; % The maximum number of packets at an SNR point

% Fixed PHY configuration for all simulations
cfgHE = wlanHEMUConfig(192); % Input 11ax allocation index 
% The full RU assignment and allocation index lookup table is shown
% in the quick start guide
% Example: when allocation index = 24, 
% then 1st RU has size 106, 2nd/3rd RU size = 52
for userIdxIter = 1:numel(cfgHE.User)
    cfgHE.User{userIdxIter}.APEPLength = 1000; % Payload length in bytes
end

%% Full PHY simulation
% Generate a structure array of simulation configurations. Each element is
% one SNR point to simulate.
simParams = getBox0SimParams(chan,numTxRx,mcs,cfgHE,maxnumberrors,maxNumPackets,...
    userIdx,ruIdx,Nsts);
snrs = [simParams.SNR];

% Simulate each configuration
results = cell(1,numel(simParams));
% Use 'parfor' to speed up the simulation
% Requirement to use 'parfor': install MATLAB Parallel Computing Toolbox 
parfor isim = 1:numel(simParams)
    results{isim} = box0Simulation(simParams(isim));
end
tEnd = toc(tStart);
fname_I = sprintf('inrPer_Config%d_%s_%s-by-%s_MCS%s.mat',cfgHE.AllocationIndex,char(chan),num2str(numTxRx(1)),num2str(numTxRx(2)),num2str(mcs));
save(fname_I,'results','mcs','numTxRx','chan','cfgHE','maxNumPackets','snrs','tEnd')

plotPERvsSNR(simParams,results);
