%% Config 192
clear all
tStart = tic;
% Simulation Parameters
mcs = [4]; % Vector of MCS to simulate between 0 and 9
numTxRx = [8 2]; % Matrix of MIMO schemes, each row is [numTx numRx]
chan = "Model-D"; % String array of delay profiles to simulate
maxnumberrors = 40*1e3;  % The maximum number of packet errors at an SNR point
maxNumPackets = 40*1e3; % The maximum number of packets at an SNR point
% maxnumberrors = 1e1;  % The maximum number of packet errors at an SNR point
% maxNumPackets = 1e2; % The maximum number of packets at an SNR point

% Fixed PHY configuration for all simulations
cfgHE = wlanHEMUConfig(97); 
for userIdx = 1:numel(cfgHE.User)
    cfgHE.User{userIdx}.APEPLength = 1000; % Payload length in bytes
end

% Generate a structure array of simulation configurations. Each element is
% one SNR point to simulate.
simParams = getBox0SimParams(chan,numTxRx,mcs,cfgHE,maxnumberrors,maxNumPackets);
snrs = [simParams.SNR];

% Simulate each configuration
results = cell(1,numel(simParams));
% parfor isim = 1:numel(simParams) % Use 'parfor' to speed up the simulation
parfor isim = 1:numel(simParams)
    results{isim} = box0Simulation1IntUser2(simParams(isim));
end
tEnd = toc(tStart);
fname_I = sprintf('snrPer_Config%d_%s_%s-by-%s_MCS%s.mat',cfgHE.AllocationIndex,char(chan),num2str(numTxRx(1)),num2str(numTxRx(2)),num2str(mcs));
save(fname_I,'results','mcs','numTxRx','chan','cfgHE','maxNumPackets','snrs','tEnd')

plotPERvsSNR(simParams,results);
