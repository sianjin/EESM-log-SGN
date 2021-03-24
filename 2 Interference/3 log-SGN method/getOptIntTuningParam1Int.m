%% Load effective SINR
clear
load('snrPer_Config192_Model-D_4-by-2_MCS4MixLoss16dB.mat')
load('eesmEffSnr_Config192_Model-D_4-by-2_MCS4MixLoss16dB.mat')
rxSnrs = snrs;
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
%% Load effective SNR
load('snrPer_Config192_Model-D_4-by-2_MCS4SigLoss16dB.mat')
load('eesmEffSnr_Config192_Model-D_4-by-2_MCS4Sig.mat')
sinrStoreSig = results{snrIdx}.sinrStore;
gammaEffdBSig = effectiveSinrVec(abstraction,sinrStoreSig,betaOpt);
gammaEffLinearSig = 10.^(gammaEffdBSig/10);
xBestSig = logSGNFitting(gammaEffLinearSig);
[effSnrSigSGNdB,perFastAbsSig] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,xBestSig);
effSnrSigSGNLinear = 10.^(effSnrSigSGNdB/10);
%% Load effective INR
load('snrPer_Config192_Model-D_4-by-2_MCS4IntLoss16dB.mat')
load('eesmEffSnr_Config192_Model-D_4-by-2_MCS4Int.mat')
sinrStoreInt = results{snrIdx}.sinrStore;
gammaEffdBInt = effectiveSinrVec(abstraction,sinrStoreInt,betaOpt);
gammaEffLinearInt = 10.^(gammaEffdBInt/10);
xBestInt = logSGNFitting(gammaEffLinearInt);
[effSnrIntSGNdB1,perFastAbsInt1] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,xBestInt);
effSnrIntSGNLinear1 = 10.^(effSnrIntSGNdB1/10); % interference 1
effSnrIntSGNLinear = {effSnrIntSGNLinear1}; % Put all interference at here
%% Model SINR distribution using log-SGN data division (not mixture model)
% thetaIni = 0.5; % Initialize interference tuning parameter 
mse = @(theta)sinrModelFittingMse(effSnrSigSGNLinear,effSnrIntSGNLinear,gammaEffdBMix,theta);
thetaOpt = fminbnd(mse,0.1,1);
% Obtain accumulative effective INR
effSnrIntSGNLinearAccum = effSnrIntSGNLinear{1};
for iter = 2:numel(effSnrIntSGNLinear)
    effSnrIntSGNLinearAccum = effSnrIntSGNLinearAccum + effSnrIntSGNLinear{iter};
end
% Obtain optimal modeled effective SINR
sinrEffModelLinear = effSnrSigSGNLinear./(1+thetaOpt*effSnrIntSGNLinearAccum);
sinrEffModeldB = 10*log10(sinrEffModelLinear);
%% EESM-log-SGN abstraction for modeled effective SINR
xBestModel = logSGNFitting(sinrEffModelLinear);
[effSinrModeldB,perFastAbsModel] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,xBestModel);
hold on
[pdfEffSnrModeldB xEffSnrModeldB]=ksdensity(effSinrModeldB);
plot(xEffSnrModeldB, pdfEffSnrModeldB)
legend('Full PHY','EESM-log-SGN-Int-Model')
xlabel('Effective SINR')
ylabel('PDF')
grid('on')
title([char(cfgHE.ChannelBandwidth) 'MHz, MCS' num2str(mcs) ', ' num2str(numTxRx(1)) 'x' num2str(numTxRx(2)) ', ' char(chan) ', ' char(coding) ', Packet Length ' num2str(dataLength) 'Byte, RX SNR ' num2str(rxSnrs(snrIdx)) 'dB']);