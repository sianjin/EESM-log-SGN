% Obtain effective SINR and average PER predicted by EESM-log-SGN
function [effSinrdB,logSGNAvgPer] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,xBest)
effSinrdB = [];
numPkts = 40000;
% Generate EESM-log-SGN effective SINR
for iter = 1: numPkts
    tempSkewGeneralizedNormRv = skewGeneralizedNormal(xBest(1), xBest(2), xBest(3), xBest(4));
    sinrEff = 10*log10(exp(tempSkewGeneralizedNormRv));
    effSinrdB = [effSinrdB, sinrEff];
end
% Obtain each PER corresponding to each effective SINR
[perVec,perPL0Vec,L0,lut] = estimatePER(abstraction,effSinrdB,format,mcs,coding,dataLength);
% Obtain averaged PER predicted by the EESM-log-SGN method
logSGNAvgPer = mean(perVec);
end