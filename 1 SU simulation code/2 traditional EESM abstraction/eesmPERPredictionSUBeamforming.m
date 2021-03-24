function out = eesmPERPredictionSUBeamforming(simParams)
% box0Simulation Example helper function

%   Copyright 2019 The MathWorks, Inc.

% Extract configuration
cfgHE = simParams.Config;
substreamidx = simParams.RandomSubstream;
maxNumPackets = simParams.MaxNumPackets;
maxNumErrors = simParams.MaxNumErrors;
snr = simParams.SNR;
beta = simParams.Beta;
% Generate channels
tgaxChannel = clone(simParams.Channel);
tgaxChannel.ChannelFiltering = false;
tgaxChannel.NormalizePathGains = true;
% Create an NDP packet with the correct number of space-time streams to
% generate enough LTF symbols
cfgNDP = wlanHESUConfig('APEPLength',0,'GuardInterval',0.8); % No data in an NDP
cfgNDP.ChannelBandwidth = cfgHE.ChannelBandwidth;
cfgNDP.NumTransmitAntennas = cfgHE.NumTransmitAntennas;
cfgNDP.NumSpaceTimeStreams = cfgHE.NumTransmitAntennas;

% Set random substream index per iteration to ensure that each
% iteration uses a repeatable set of random numbers
stream = RandStream('combRecursive','Seed',99);
stream.Substream = substreamidx;
RandStream.setGlobalStream(stream);


% Indices to extract fields from the PPDU
ind = wlanFieldIndices(cfgHE);

% Get occupied subcarrier indices and OFDM parameters
ofdmInfo = wlanHEOFDMInfo('HE-Data',cfgHE);

% Create an instance of the AWGN channel per SNR point simulated
awgnChannel = comm.AWGNChannel;
awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
% Account for noise energy in nulls so the SNR is defined per
% active subcarrier
awgnChannel.SNR = snr-10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);
N0 = 10^(-awgnChannel.SNR/10); 

% Get path filers for last channel (same for all channels)
chanInfo = info(tgaxChannel);
pathFilters = chanInfo.ChannelFilterCoefficients; % [NP numChannelTaps]

% Create object to deal with abstraction
Abstraction = tgaxEESMLinkPerformanceModel;

% Loop to simulate multiple packets
numPacketErrors = 0;
numPacketErrorsAbs = 0;
numPkt = 1; % Index of packet transmitted
% Storing effective SNR vector for estimating effective SNR distribution
effSnrVec = [];
while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
    reset(tgaxChannel); % Reset channel for different realization
    
    pathGains = tgaxChannel();    % Get path gains
    chan = helperPerfectChannelEstimate(pathGains,pathFilters,ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices);
    
    % -------------------------------------------------------------------
    % For each user STA, pass the NDP packet through the channel and calculate
    % the feedback channel state matrix by SVD.
    staFeedback = heIdealUserBeamformingFeedback(tgaxChannel,cfgNDP); % Get the full-band beamforming feedback for a user
    steeringMatrix = heSUCalculateSteeringMatrix(staFeedback,cfgHE,cfgNDP);
    cfgHE.SpatialMapping = 'Custom';
    cfgHE.SpatialMappingMatrix = steeringMatrix;
    % -------------------------------------------------------------------
    
    % Per user processing, focus on the first user
    Nsts = cfgHE.NumSpaceTimeStreams;
    % Calculate SINR using abstraction
    % As multiple symbols returned average over symbols and permute
    % for calculations
    % Get precoding matrix for abstraction
    Wtx = getPrecodingMatrix(cfgHE); % Include cyclic shift applied per STS
    WtxUser= Wtx/sqrt(Nsts);
    Htxrx = permute(mean(chan,2),[1 3 4 2]); % Nst-by-Nt-by-Nr
    Ptxrx = 1; % Assume transmit power is 0dBW
    sinr = calculateSINR(Htxrx,Ptxrx,WtxUser,N0);
    
    % Link performance model - estimate PER using abstraction
    snrEff = eesmEffectiveSINR(Abstraction,sinr,beta);
    effSnrVec = [effSnrVec, snrEff];
    perIns = estimateLinkPerformance(Abstraction,snrEff,cfgHE);

    % Flip a coin for the abstracted PHY
    packetErrorAbs = rand(1)<=perIns;
    numPacketErrorsAbs = numPacketErrorsAbs+packetErrorAbs;

    numPkt = numPkt+1;
end

% Remove last increment
numPkt = numPkt-1;

% Calculate packet error rate (PER) at SNR point
packetErrorRateAbs = numPacketErrorsAbs/numPkt;

% Return results
out = struct;
out.packetErrorRateAbs = packetErrorRateAbs;
out.effSnrVec = effSnrVec;

disp([char(simParams.DelayProfile) ' '...
      num2str(simParams.NumTransmitAntennas) '-by-' ...
      num2str(simParams.NumReceiveAntennas) ','...
      ' MCS ' num2str(simParams.MCS) ','...
      ' SNR ' num2str(simParams.SNR) ...
      ' completed after ' num2str(numPkt) ' packets,'...
      ' PER:' num2str(packetErrorRateAbs)]);
  
end