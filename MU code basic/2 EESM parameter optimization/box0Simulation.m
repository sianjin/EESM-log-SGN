function effSnrVec = box0Simulation(simParams)
% box0Simulation Example helper function

%   Copyright 2019 The MathWorks, Inc.

% Per user processing
userIdx = 1;

% Extract configuration
cfgHE = simParams.Config;
substreamidx = simParams.RandomSubstream;
maxNumPackets = simParams.MaxNumPackets;
maxNumErrors = simParams.MaxNumErrors;
snr = simParams.SNR;
tgaxChannel = simParams.Channel;
tgaxChannel.UserIndex = userIdx;



% Set random substream index per iteration to ensure that each
% iteration uses a repeatable set of random numbers
stream = RandStream('combRecursive','Seed',99);
stream.Substream = substreamidx;
RandStream.setGlobalStream(stream);

% Indices to extract fields from the PPDU
ind = wlanFieldIndices(cfgHE);

% Get occupied subcarrier indices and OFDM parameters
ofdmInfo = wlanHEOFDMInfo('HE-Data',cfgHE,userIdx);

% Create an instance of the AWGN channel per SNR point simulated
awgnChannel = comm.AWGNChannel;
awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
% Account for noise energy in nulls so the SNR is defined per
% active subcarrier
awgnChannel.SNR = snr-10*log10(ofdmInfo.FFTLength/sum(ruInfo(cfgHE).RUSizes));

% For abstraction
Wtx = getPrecodingMatrix(cfgHE,userIdx); % Include cyclic shift applied per STS
Wtx = Wtx/sqrt(cfgHE.User{userIdx}.NumSpaceTimeStreams);
% N0 = 10^(-snr/10); original but I don't think it's correct
N0 = 10^(-awgnChannel.SNR/10); 

% Get path filers for last channel (same for all channels)
chanInfo = info(tgaxChannel);
pathFilters = chanInfo.ChannelFilterCoefficients; % [NP numChannelTaps]

% Create object to deal with abstraction
Abstraction = tgaxLinkPerformanceModel;

% Loop to simulate multiple packets
perStore = nan(maxNumPackets,1);
perAbsStore = nan(maxNumPackets,1);
perAbsRawStore = nan(maxNumPackets,1);
snreffStore = nan(maxNumPackets,1);
numPacketErrors = 0;
numPacketErrorsAbs = 0;
numPkt = 1; % Index of packet transmitted

% storing effective SNR vector for estimating effective SNR distribution
effSnrVec = [];

while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
    % Generate a packet with random PSDU
    psduLengthVec = getPSDULength(cfgHE); % PSDU length in bytes
    psduLength = psduLengthVec(userIdx); % hard code to the psduLength of the 1st user
    txPSDU = randi([0 1],psduLength*8,1,'int8');
    tx = wlanWaveformGenerator(txPSDU,cfgHE);
    
    % Add trailing zeros to allow for channel delay
    txPad = [tx; zeros(50,cfgHE.NumTransmitAntennas)];

    % Pass through a fading indoor TGax channel
    reset(tgaxChannel); % Reset channel for different realization
    [rx,pathGains] = tgaxChannel(txPad);

    % Get perfect timing offset and channel matrix for HE-LTF field
    heltfPathGains = pathGains(ind.HELTF(1):ind.HELTF(2),:,:,:,:);
    pktOffset = channelDelay(heltfPathGains,pathFilters);
    chan = helperPerfectChannelEstimate(heltfPathGains,pathFilters,ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices,pktOffset);

    % Calculate SINR using abstraction
    % As multiple symbols returned average over symbols and permute
    % for calculations
    Htxrx = permute(mean(chan,2),[1 3 4 2]); % Nst-by-Nt-by-Nr
    Ptxrx = 1; % Assume transmit power is 0dBW
    sinr = calculateSINR(Htxrx,Ptxrx,Wtx,N0);

    % Link performance model - estimate PER using abstraction
    [perAbs,effSINR] = estimateLinkPerformance(Abstraction,sinr,cfgHE);
    
    effSnrVec = [effSnrVec, effSINR];

    numPkt = numPkt+1;
end
  
end