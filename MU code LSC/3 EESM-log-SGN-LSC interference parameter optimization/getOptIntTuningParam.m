%% Introduction 
% This file is the main script to optimize EESM-log-SGN-LSC
% interference parameter theta that is used in equation (6) 
% of the IEEE TCOM paper:
% "Efficient PHY Layer Abstraction for Fast Simulations in Complex 
% System Environments"
% * If the simulation stops with a bug, re-run this file would be fine.
%% Load basic setup
clear all
coding = 'LDPC';
dataLength = 1000;
format = 'HE_MU';
abstraction = tgaxEESMLinkPerformanceModel;
%% Load effective SNR
load('snrPer_Config192_Model-D_4-by-2_MCS4.mat')
load('eesmEffSnr_Config192_Model-D_4-by-2_MCS4.mat')
snrIdxSig = 8;
snrStore = results{snrIdxSig}.sinrStore;
effSNRdB = effectiveSinrVec(abstraction,snrStore,betaOptSig);
effSNRLinear = 10.^(effSNRdB/10);
logSGNParamSig = logSGNFitting(effSNRLinear);
[effSNRSGNdB,logSGNAvgPerSig] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,logSGNParamSig);
effSNRSGNLinear = 10.^(effSNRSGNdB/10);
%% Load effective INR for the 1st interferer
load('inrPer_Config192_Model-D_4-by-2_MCS4.mat')
load('eesmEffSnr_Config192_Model-D_4-by-2_MCS4.mat')
inrIdx = 4;
inrStore = results{inrIdx}.sinrStore;
effINRdB1 = effectiveSinrVec(abstraction,inrStore,betaOptSig);
effINRLinear1 = 10.^(effINRdB1/10);
logSGNParamInt1 = logSGNFitting(effINRLinear1);
[effINRSGNdB1,logSGNAvgPerInt1] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,logSGNParamInt1);
effINRSGNLinear1 = 10.^(effINRSGNdB1/10);  
%% Form effective INR cell if there are multiple interferers
effINRSGNLinear = {effINRSGNLinear1}; % If there are multiple interference, put them in this cell
%% Load effective SINR for calibration
load('sinrPer_Config192_Model-D_4-by-2_MCS4.mat')
load('eesmEffSinr_Config192_Model-D_4-by-2_MCS4.mat')
rxSnrs = snrs;
snrIdx = 4;
sinrStore = results{snrIdx}.sinrStore;
abstraction = tgaxEESMLinkPerformanceModel;
effSINRdB = effectiveSinrVec(abstraction,sinrStore,betaOpt);
%% True effective SINR histogram
histogram(effSINRdB, 'normalization', 'pdf')
%% Optimize interference parameter theta used in equation (6) of the IEEE TCOM paper
mse = @(theta)sinrModelFittingMse(effSNRSGNLinear,effINRSGNLinear,effSINRdB,theta);
thetaOpt = fminbnd(mse,0.1,1);
% Obtain accumulative effective INR in the denomenator of in equation (6)
effSnrIntSGNLinearAccum = effINRSGNLinear{1};
for iter = 2:numel(effINRSGNLinear)
    effSnrIntSGNLinearAccum = effSnrIntSGNLinearAccum + effINRSGNLinear{iter};
end
% Calculate equation (6)
effSINRLSCLinear = effSNRSGNLinear./(1+thetaOpt*effSnrIntSGNLinearAccum);
%% EESM-log-SGN abstraction for modeled effective SINR
logSGNParamLSC = logSGNFitting(effSINRLSCLinear);
[effSINRLSCdB,logSGNAvgPerLSC] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,logSGNParamLSC);
hold on
[pdfEffSnrModeldB xEffSnrModeldB]=ksdensity(effSINRLSCdB);
plot(xEffSnrModeldB, pdfEffSnrModeldB,'LineWidth',2)
legend('Traditional EESM','EESM-log-SGN-LSC')
xlabel('Effective SINR')
ylabel('PDF')
grid('on')
title(['MCS' num2str(mcs) ', ' num2str(numTxRx(1)) 'x' num2str(numTxRx(2)) ', ' char(chan) ', ' char(coding) ', Packet Length ' num2str(dataLength) 'Byte, RX SNR ' num2str(rxSnrs(snrIdx)) 'dB']);