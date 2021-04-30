%% Introduction 
% This file is the main script to run traditional EESM PHY layer
% abstraction shown in Fig.1 in the IEEE TCOM paper:
% "Efficient PHY Layer Abstraction for Fast Simulations in Complex 
% System Environments"
% The traditional EESM simulation configurations in this script cover: 
% 1) 20MHz OFDM/OFDMA MIMO/MU-MIMO system
% 2) No interference 
%% Traditional EESM PHY layer abstraction setup
clear all
tStart = tic;
mcs = [4]; % Vector of MCS to simulate between 0 and 11
numTxRx = [8 2]; % Matrix of MIMO schemes, each row is [numTx numRx]
chan = "Model-D"; % String array of delay profiles to simulate
userIdx = 2; % User of investigation
ruIdx = 2; % RU of investigation
Nsts = 2; % Number of space-time streams
betaOpt = 8.3891; % EESM tuning parameter
maxnumberrors = 40*1e3;  % The maximum number of packet errors at an SNR point
maxNumPackets = 40*1e3; % The maximum number of packets at an SNR point
% maxnumberrors = 5*1e1;  % The maximum number of packet errors at an SNR point
% maxNumPackets = 5*1e2; % The maximum number of packets at an SNR point

% Fixed PHY configuration for all simulations
cfgHE = wlanHEMUConfig(97); 
for userIdxIter = 1:numel(cfgHE.User)
    cfgHE.User{userIdxIter}.APEPLength = 1000; % Payload length in bytes
end

%% Traditional EESM PHY layer abstraction
% Generate a structure array of simulation configurations. Each element is
% one SNR point to simulate.
simParams = getBox0SimParams(chan,numTxRx,mcs,cfgHE,maxnumberrors,maxNumPackets,...
    userIdx,ruIdx,Nsts,betaOpt);
snrs = [simParams.SNR];

% Simulate each configuration
avgPerEESM = zeros(1,numel(simParams));
results = cell(1,numel(simParams));
% parfor isim = 1:numel(simParams) % Use 'parfor' to speed up the simulation
parfor isim = 1:numel(simParams)
    out = box0Simulation(simParams(isim));
    results{isim} = out;
    avgPerEESM(isim) = out.packetErrorRateAbs;
end
tEndEesmAbs = toc(tStart);
fname_I = sprintf('eesmAvgPer_Config%d_%s_%s-by-%s_MCS%s.mat',cfgHE.AllocationIndex,char(chan),num2str(numTxRx(1)),num2str(numTxRx(2)),num2str(mcs));
save(fname_I,'results','mcs','numTxRx','chan','cfgHE','maxNumPackets','snrs','tEndEesmAbs')
plot(snrs, avgPerEESM,'-o','LineWidth',1)
title('Traditional EESM PHY Layer Abstraction')
xlabel('RX SNR (dB)');
ylabel('Average PER');
set(gca, 'YScale', 'log')
legend('EESM, MCS4')
grid on