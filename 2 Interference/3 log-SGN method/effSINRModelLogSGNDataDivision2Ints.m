% Simulation might stop and pop up with an error. 
% This is caused by inproper initial values in functions.
% Rerun this code can solve this issue.
%% Load effective SINR
clear
load('snrPer_Config192_Model-D_8-by-2_MCS4MixLoss12dBInt2.mat')
rxSnrs = snrs;
load('eesmEffSnr_Config192_Model-D_8-by-2_MCS4MixLoss12dBInt2.mat')
coding = 'LDPC';
dataLength = 1000;
format = 'HE_MU';
snrIdx = 4;
sinrStore = results{snrIdx}.sinrStore;
abstraction = tgaxEESMLinkPerformanceModel;
gammaEffdBMix = effectiveSinrVec(abstraction,sinrStore,betaOpt);
gammaEffLinearMix = 10.^(gammaEffdBMix/10);
%% True SINR histogram
histogram(gammaEffdBMix, 'normalization', 'pdf')
%% EESM-log-SGN abstraction for true effective SINR
xBestMix = logSGNFitting(gammaEffLinearMix);
tStartSINR = tic;
[effSinrdB,perFastAbs] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,xBestMix);
tEndSINR = toc(tStartSINR);
hold on
[pdfEffSnrdB xEffSnrdB]=ksdensity(effSinrdB);
plot(xEffSnrdB, pdfEffSnrdB)
%% Load effective SNR
load('snrPer_Config192_Model-D_8-by-2_MCS4SigLoss12dB.mat')
load('eesmEffSnr_Config192_Model-D_8-by-2_MCS4Sig.mat')
sinrStoreSig = results{snrIdx}.sinrStore;
gammaEffdBSig = effectiveSinrVec(abstraction,sinrStoreSig,betaOpt);
gammaEffLinearSig = 10.^(gammaEffdBSig/10);
xBestSig = logSGNFitting(gammaEffLinearSig);
%% Load effective INR
load('snrPer_Config192_Model-D_8-by-2_MCS42IntsLoss12dB.mat')
load('eesmEffSnr_Config192_Model-D_8-by-2_MCS4Int.mat')
sinrStoreInt = results{snrIdx}.sinrStore;
gammaEffdBInt = effectiveSinrVec(abstraction,sinrStoreInt,betaOpt);
gammaEffLinearInt = 10.^(gammaEffdBInt/10);
xBestInt = logSGNFitting(gammaEffLinearInt);
%% Model SINR distribution using log-SGN data division (not mixture model)
tStartSINRModel = tic;
[effSnrSigSGNdB,perFastAbsSig] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,xBestSig);
effSnrSigSGNLinear = 10.^(effSnrSigSGNdB/10);
[effSnrIntSGNdB,perFastAbsInt] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,xBestInt);
effSnrIntSGNLinear = 10.^(effSnrIntSGNdB/10);
[effSnrInt2SGNdB,perFastAbsInt2] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,xBestInt);
effSnrInt2SGNLinear = 10.^(effSnrInt2SGNdB/10);
theta = 0.7206;
sinrEffModelLinear = effSnrSigSGNLinear./(1+theta*(effSnrIntSGNLinear+effSnrInt2SGNLinear));
sinrEffModeldB = 10*log10(sinrEffModelLinear);
tEndSINRModel = toc(tStartSINRModel);
% hold on
% histogram(sinrEffModeldB, 'normalization', 'pdf')
%% EESM-log-SGN abstraction for modeled effective SINR
xBestModel = logSGNFitting(sinrEffModelLinear);
[effSinrModeldB,perFastAbsModel] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,xBestModel);
hold on
[pdfEffSnrModeldB xEffSnrModeldB]=ksdensity(effSinrModeldB);
plot(xEffSnrModeldB, pdfEffSnrModeldB)
legend('Full PHY','EESM-log-SGN', 'EESM-log-SGN-Int-Model')
xlabel('Effective SINR')
ylabel('PDF')
grid('on')
title([char(cfgHE.ChannelBandwidth) 'MHz, MCS' num2str(mcs) ', ' num2str(numTxRx(1)) 'x' num2str(numTxRx(2)) ', ' char(chan) ', ' char(coding) ', Packet Length ' num2str(dataLength) 'Byte, RX SNR ' num2str(rxSnrs(snrIdx)) 'dB']);