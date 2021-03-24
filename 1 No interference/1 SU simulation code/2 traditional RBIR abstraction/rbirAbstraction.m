clear all
tStart = tic;
% Simulation Parameters
mcs = [4]; % Vector of MCS to simulate between 0 and 9
numTxRx = [1 1]; % Matrix of MIMO schemes, each row is [numTx numRx]
chan = "Model-D"; % String array of delay profiles to simulate
maxnumberrors = 40*1e3;  % The maximum number of packet errors at an SNR point
maxNumPackets = 40*1e3; % The maximum number of packets at an SNR point
% maxnumberrors = 1e1;  % The maximum number of packet errors at an SNR point
% maxNumPackets = 1e2; % The maximum number of packets at an SNR point
rbirParameters = [1.0355, 1.0726]; % RBIR tuning parameter

% Fixed PHY configuration for all simulations
cfgHE = wlanHEMUConfig(192); 
for userIdx = 1:numel(cfgHE.User)
    cfgHE.User{userIdx}.APEPLength = 1000; % Payload length in bytes
end

% Generate a structure array of simulation configurations. Each element is
% one SNR point to simulate.
simParams = getBox0SimParams(chan,numTxRx,mcs,cfgHE,maxnumberrors,maxNumPackets,rbirParameters);
snrs = [simParams.SNR];

% Simulate each configuration
avgPerRbir = zeros(1,numel(simParams));
results = cell(1,numel(simParams));
% parfor isim = 1:numel(simParams) % Use 'parfor' to speed up the simulation
parfor isim = 1:numel(simParams)
    out = rbirPERPredictionMUBeamforming(simParams(isim));
    results{isim} = out;
    avgPerRbir(isim) = out.packetErrorRateAbs;
end
tEndRbirAbs = toc(tStart);
fname_I = sprintf('rbirTxSnrAvgPer_Config%d_%s_%s-by-%s_MCS%s.mat',cfgHE.AllocationIndex,char(chan),num2str(numTxRx(1)),num2str(numTxRx(2)),num2str(mcs));
save(fname_I,'results','mcs','numTxRx','chan','cfgHE','maxNumPackets','snrs','avgPerRbir','tEndRbirAbs')
plot(snrs, avgPerRbir)

