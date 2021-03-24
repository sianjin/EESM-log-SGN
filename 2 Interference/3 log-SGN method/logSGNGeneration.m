function [effSinrdB,perFastAbs] = logSGNGeneration(abstraction,coding,mcs,dataLength,format,xBest)
effSinrdB = [];
numPacketErrorsAbs = 0;
numPkts = 40000;
for iter = 1: numPkts
    tempSkewGeneralizedNormRv = skewGeneralizedNormal(xBest(1), xBest(2), xBest(3), xBest(4));
    sinrEff = 10*log10(exp(tempSkewGeneralizedNormRv));
    effSinrdB = [effSinrdB, sinrEff];
%     perIns = estimatePER(abstraction,sinrEff,format,mcs,coding,dataLength);
%     packetErrorAbs = rand(1)<=perIns;
%     numPacketErrorsAbs = numPacketErrorsAbs+packetErrorAbs;
end
[perVec,perPL0Vec,L0,lut] = estimatePER(abstraction,effSinrdB,format,mcs,coding,dataLength);
perFastAbs = mean(perVec);
% perFastAbs = numPacketErrorsAbs/numPkts;
end