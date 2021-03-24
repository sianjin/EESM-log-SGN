% calculate autocorrelation function of the effective SNR 
% when channel is correlated and generated using 
% eesmPERPredictionSUBeamformingCorrelation.m
load('eesmTxSnrAvgPer_CBW20_Model-D_1-by-1_MCS4IID.mat');
effSnrProcessIID = results{1}.effSnrVec;
meanEffSnr = mean(effSnrProcessIID); % mean of the random process
varEffSnr = var(effSnrProcessIID); % var of the random process
load('eesmTxSnrAvgPer_CBW20_Model-D_1-by-1_MCS4.mat');
effSnrProcess = results{1}.effSnrVec;
n = 40000; % length of the random process
autoCorrEffSnr = []; 
for lag = 1:100 % 1 to max lag 
    autoCorrEffSnr_lag = 0; % autocorrelation function at each lag    
    for t = 1:n-lag
        autoCorrEffSnr_lag = autoCorrEffSnr_lag + (effSnrProcess(t) - meanEffSnr)*(effSnrProcess(t+lag) - meanEffSnr);
    end
    autoCorrEffSnr_lag = autoCorrEffSnr_lag/(n-lag)/varEffSnr;
    autoCorrEffSnr = [autoCorrEffSnr, autoCorrEffSnr_lag];
end
plot(autoCorrEffSnr)