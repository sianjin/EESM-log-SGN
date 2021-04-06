%% Introduction
% This file is the main script to run the EESM-log-SGN mixture model
% corresponding to Algorithm 3 in the IEEE TCOM paper:
% "Efficient PHY Layer Abstraction for Fast Simulations in Complex 
% System Environments"
%% Load setup
coding = 'LDPC'; % Channel coding
dataLength = 1000; % APEP length 
format = 'HE_MU'; % OFDM/OFDMA MIMO/MU-MIMO setup
mcs = 4; % MCS 
eesmAbstraction = tgaxEESMLinkPerformanceModel;
snr1 = 20; % RX SNR 1 in dB
snr2 = 22; % RX SNR 2 in dB
snr = 21; % RX SNR under investigation in dB
% Optimized log-SGN parameters under RX SNR 1
logSGNParamBest1 = [
2.5712    0.3219    1.7736    0.0000
];
% Optimized log-SGN parameters under RX SNR 2
logSGNParamBest2 = [
2.8356    0.3418    1.5308    0.0000
];
%% Algorithm 3 in our IEEE TCOM paper
effSnrdB = []; % Effective SNR histogram the EESM-log-SGN mixture model
numPkts = 40000; % Number of packets
for iter = 1: numPkts
    epsilon = (snr-snr1)/(snr2-snr1); % Mixture probability
    u = rand;
    if u < 1 - epsilon
        % Generate log-SGN random variable under RX SNR 1
        tempSkewGeneralizedNormRv1 = skewGeneralizedNormal(logSGNParamBest1(1), logSGNParamBest1(2), logSGNParamBest1(3), logSGNParamBest1(4));
        snrEff1 = 10*log10(exp(tempSkewGeneralizedNormRv1));
        snrEff = snrEff1; % Effective SNR in dB
    else
        % Generate log-SGN random variable under RX SNR 2
        tempSkewGeneralizedNormRv2 = skewGeneralizedNormal(logSGNParamBest2(1), logSGNParamBest2(2), logSGNParamBest2(3), logSGNParamBest2(4));
        snrEff2 = 10*log10(exp(tempSkewGeneralizedNormRv2));
        snrEff = snrEff2; % Effective SNR in dB
    end
    effSnrdB = [effSnrdB, snrEff]; 
end
%% Plot effective SNR histogram estimated by the EESM-log-SGN mixture model
[pdfEffSnrdB xEffSnrdB]=ksdensity(effSnrdB);
plot(xEffSnrdB, pdfEffSnrdB,'--','LineWidth',2)
grid('on')
xlabel('Effective SINR (dB)')
ylabel('PDF')
legend('EESM-log-SGN Mixture')
%% Obtain average PER estimated from the EESM-log-SGN mixture model
[perVec,perPL0Vec,L0,lut] = estimatePER(eesmAbstraction,effSnrdB,format,mcs,coding,dataLength);
% Average PER estimated from the EESM-log-SGN mixture model
logSGNMixtureAvgPer = mean(perVec);

