cfgHE = wlanHEMUConfig(5);
cfgHE.SIGBCompression = false;
cfgHE.NumTransmitAntennas = 1;

psdu = cell(1,numel(cfgHE.User));
psduLength = getPSDULength(cfgHE);
for j = 1:numel(cfgHE.User)
    psdu = randi([0 1],psduLength(j)*8,1,'int8');
end

y = wlanWaveformGenerator(psdu,cfgHE);
plot(abs(y))

% cfgHE = wlanHESUConfig;
% cfgHE.NumTransmitAntennas = 1;
% 
% psduLength = getPSDULength(cfgHE);
% psdu = randi([0 1],psduLength*8,1,'int8');
% 
% y = wlanWaveformGenerator(psdu,cfgHE);
% plot(abs(y))