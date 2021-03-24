effSinrdB = [];
numPacketErrorsAbs = 0;
numPkts = 40000;
coding = 'LDPC';
dataLength = 1000;
format = 'HE_MU';
mcs = 4;
eesmAbstraction = tgaxEESMLinkPerformanceModel;
theta = 0.4421;
xBestSig1 = [
 2.3875    0.1771   -1.0050    2.9028
];
xBestInt1 = [
 -3.0184    0.4331   -0.2396    2.0526
];
xBestSig2 = [
 2.6060    0.2088   -0.8122    2.3540
];
xBestInt2 = [
 -2.5541    0.4340   -0.2551    1.6913
];
snr1 = 16;
snr2 = 18;
snr = 17;
for iter = 1: numPkts
    p = (snr-snr1)/(snr2-snr1); % mixture probability
    if binornd(1,p) == 0
        tempSGNSig = skewGeneralizedNormal(xBestSig1(1), xBestSig1(2), xBestSig1(3), xBestSig1(4));
        tempSGNInt = skewGeneralizedNormal(xBestInt1(1), xBestInt1(2), xBestInt1(3), xBestInt1(4));
    else
        tempSGNSig = skewGeneralizedNormal(xBestSig2(1), xBestSig2(2), xBestSig2(3), xBestSig2(4));
        tempSGNInt = skewGeneralizedNormal(xBestInt2(1), xBestInt2(2), xBestInt2(3), xBestInt2(4));
    end
    snrEffSigLinear = exp(tempSGNSig);
    snrEffIntLinear = exp(tempSGNInt);
    snrEff = snrEffSigLinear/(1+theta*snrEffIntLinear);
    snrEff = 10*log10(snrEff);
    effSinrdB = [effSinrdB, snrEff];
end
[perVec,perPL0Vec,L0,lut] = estimatePER(eesmAbstraction,effSinrdB,format,mcs,coding,dataLength);
perFastAbs = mean(perVec);
hold on
[pdfEffSnrdB xEffSnrdB]=ksdensity(effSinrdB);
plot(xEffSnrdB, pdfEffSnrdB)
