%% get parameters
load('snrPer_Config24_Model-D_8-by-2_MCS4MixLoss20dB.mat')
load('eesmEffSnr_Config24_Model-D_8-by-2_MCS4IntLoss20dB.mat')
coding = 'LDPC';
dataLength = 1000;
format = 'HE_MU';
snrIdx = 2;
sinrStore = results{snrIdx}.sinrStore;
abstraction = tgaxEESMLinkPerformanceModel;
gammaEffdB = effectiveSinrVec(abstraction,sinrStore,betaOpt);
gammaEffLinear = 10.^(gammaEffdB/10);
logGammaEffLinear = log(gammaEffLinear);
n = length(logGammaEffLinear);
M = mean(logGammaEffLinear);
S = std(logGammaEffLinear);
xBest = zeros(1,4);
fBest = -10e10;
randRange1 = 1;
randRange2 = 0.1;
x0 = [M,S, -randRange1 + (2*randRange1)*rand, (2*randRange2)*rand];
A = [0, 0, 0, -1; 0, -1, 0, 0];
b = [0, 0];
for iter = 1: 100
    LikelihoodFun = @(x) -(-n/2* log(pi * x(2)^2/2) - 1/2/x(2)^2*sum((logGammaEffLinear-x(1)).^2) + sum(log(normcdf(x(3).*(logGammaEffLinear-x(1))./sqrt(x(2)^2+ x(4).*(logGammaEffLinear-x(1)).^2)))));
    [x, fval] = fmincon(LikelihoodFun,x0,A,b);
    if fval>fBest
        fBest = fval;
        xBest = x;
        x0 = xBest;
    end
end
%% generate and plot approximated effective SNR distribution
effSnrdB = [];
numPacketErrorsAbs = 0;
numPkts = 40000;
tic;
for iter = 1: numPkts
    tempSkewGeneralizedNormRv = skewGeneralizedNormal(xBest(1), xBest(2), xBest(3), xBest(4));
    snrEff = 10*log10(exp(tempSkewGeneralizedNormRv));
    effSnrdB = [effSnrdB, snrEff];
    perIns = estimatePER(abstraction,snrEff,format,mcs,coding,dataLength);
    packetErrorAbs = rand(1)<=perIns;
    numPacketErrorsAbs = numPacketErrorsAbs+packetErrorAbs;
end
perFastAbs = numPacketErrorsAbs/numPkts;
tEndAbs = toc;
histogram(gammaEffdB, 'normalization', 'pdf')
hold on
[pdfEffSnrdB xEffSnrdB]=ksdensity(effSnrdB);
plot(xEffSnrdB, pdfEffSnrdB)
xlabel('Effective SINR')
ylabel('PDF')
legend('Full PHY','EESM-log-SGN')
title([char(cfgHE.ChannelBandwidth) 'MHz, MCS' num2str(mcs) ', ' num2str(numTxRx(1)) 'x' num2str(numTxRx(2)) ', ' char(chan) ', ' char(coding) ', Packet Length ' num2str(dataLength) 'Byte, RX SNR ' num2str(snrs(snrIdx)) 'dB']);
grid('on')