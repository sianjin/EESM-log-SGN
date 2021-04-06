%% Introduction 
% This file is the main script to optimize EESM-log-SGN-LSC
% in Fig. 11 of the IEEE TCOM paper:
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
snrIdxSig = 3;
snrStore = results{snrIdxSig}.sinrStore;
effSNRdB = effectiveSinrVec(abstraction,snrStore,betaOptSig);
effSNRLinear = 10.^(effSNRdB/10);
logSGNParamSig = logSGNFitting(effSNRLinear);
%% Load effective INR for the 1st interferer
load('inrPer_Config192_Model-D_4-by-2_MCS4.mat')
load('eesmEffSnr_Config192_Model-D_4-by-2_MCS4.mat')
inrIdx = 1;
inrStore = results{inrIdx}.sinrStore;
effINRdB1 = effectiveSinrVec(abstraction,inrStore,betaOptSig);
effINRLinear1 = 10.^(effINRdB1/10);
logSGNParamInt1 = logSGNFitting(effINRLinear1); 
%% Load effective SINR for calibration
load('sinrPer_Config192_Model-D_4-by-2_MCS4.mat')
load('eesmEffSinr_Config192_Model-D_4-by-2_MCS4.mat')
rxSnrs = snrs;
snrIdx = 1;
sinrStore = results{snrIdx}.sinrStore;
abstraction = tgaxEESMLinkPerformanceModel;
effSINRdB = effectiveSinrVec(abstraction,sinrStore,betaOpt);
%% True effective SINR histogram
histogram(effSINRdB, 'normalization', 'pdf')
%% Model SINR distribution using EESM-log-SGN-LSC
tStartLSC = tic;
[effSNRSGNdB,logSGNAvgPerSig] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,logSGNParamSig);
effSNRSGNLinear = 10.^(effSNRSGNdB/10);
[effINRSGNdB1,logSGNAvgPerInt1] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,logSGNParamInt1);
effINRSGNLinear1 = 10.^(effINRSGNdB1/10); 
thetaOpt = 0.6434;
effSINRLSCLinear = effSNRSGNLinear./(1+thetaOpt*(effINRSGNLinear1));
tEndLSC = toc(tStartLSC);
%% EESM-log-SGN abstraction for modeled effective SINR
logSGNParamLSC = logSGNFitting(effSINRLSCLinear);
[effSINRSGNLSCdB,logSGNAvgPerLSC] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,logSGNParamLSC);
hold on
[pdfeffSINRSGNLSCdB xeffSINRSGNLSCdB]=ksdensity(effSINRSGNLSCdB);
plot(xeffSINRSGNLSCdB, pdfeffSINRSGNLSCdB,'LineWidth',2)
legend('Traditional EESM','EESM-log-SGN-LSC')
xlabel('Effective SINR')
ylabel('PDF')
grid('on')
title(['MCS' num2str(mcs) ', ' num2str(numTxRx(1)) 'x' num2str(numTxRx(2)) ', ' char(chan) ', ' char(coding) ', Packet Length ' num2str(dataLength) 'Byte, RX SNR ' num2str(rxSnrs(snrIdx)) 'dB']);