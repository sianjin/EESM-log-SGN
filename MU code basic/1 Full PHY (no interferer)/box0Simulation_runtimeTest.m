% Interference-free full PHY simulation
% under 20MHz OFDM/OFDMA MIMO/MU-MIMO system
function out = box0Simulation_runtimeTest(simParams)
% box0Simulation Example helper function

% Copyright 2019 The MathWorks, Inc.

% Extract configuration
userIdx = simParams.userIdx;
ruIdx = simParams.ruIdx;
cfgHE = simParams.Config;
substreamidx = simParams.RandomSubstream;
maxNumPackets = simParams.MaxNumPackets;
maxNumErrors = simParams.MaxNumErrors;
numUsers = size(cfgHE.User,2);
snr = simParams.SNR;
% Generate per-user channels, adapted from HEDownlinkMUExample
tgaxChannel = cell(1,numUsers);
% Generate per-user channels
for userIdxIter = 1:numUsers
    tgaxChannel{userIdxIter} = clone(simParams.Channel);
    tgaxChannel{userIdxIter}.UserIndex = userIdxIter; % Set unique user index
end
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

% Number of space-time streams
Nsts = cfgHE.User{userIdx}.NumSpaceTimeStreams;
% Get occupied subcarrier indices and OFDM parameters
ofdmInfo = wlanHEOFDMInfo('HE-Data',cfgHE,ruIdx);

% Create an instance of the AWGN channel per SNR point simulated
awgnChannel = comm.AWGNChannel;
awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
% Account for noise energy in nulls so the SNR is defined per
% active subcarrier
awgnChannel.SNR = snr-10*log10(ofdmInfo.FFTLength/sum(ruInfo(cfgHE).RUSizes));
N0 = 10^(-awgnChannel.SNR/10); 

% Get path filers for last channel (same for all channels)
chanInfo = info(tgaxChannel{userIdx});
pathFilters = chanInfo.ChannelFilterCoefficients; % [NP numChannelTaps]

% Loop to simulate multiple packets
perStore = nan(maxNumPackets,1);
numPacketErrors = 0;
numPkt = 1; % Index of packet transmitted
while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
    % Reset channel for different realization
    for userIdxIter = 1:numel(tgaxChannel)
        reset(tgaxChannel{userIdxIter}); 
    end
    
    % The below dashed part caluclates ZF precoding matrix  
    % Comment the below dashed part to run simulation with default Fourier precoding
    % -------------------------------------------------------------------
    % Generate NDP packet - with an empty PSDU as no data
    txNDP = wlanWaveformGenerator([],cfgNDP);
    % For each user STA, pass the NDP packet through the channel and calculate
    % the feedback channel state matrix by SVD
    staFeedback = cell(1,numUsers);
    for userIdxIter = 1:numel(tgaxChannel)
        % Received waveform at user STA with 50 sample padding. No noise.
        rxNDP = tgaxChannel{userIdxIter}([txNDP; zeros(50,size(txNDP,2))]);
        
        % Get the full-band beamforming feedback for a user
        staFeedback{userIdxIter} = heUserBeamformingFeedback(rxNDP,cfgNDP);
    end
    % Calculate the steering matrix to apply to the RU given the feedback
    % A zero forcing solution is used to calculate the steering matrix
    steeringMatrix = heMUCalculateSteeringMatrix(staFeedback,cfgHE,cfgNDP,ruIdx);
    
    % Apply the steering matrix to the RU of interest
    cfgHE.RU{ruIdx}.SpatialMapping = 'Custom';
    cfgHE.RU{ruIdx}.SpatialMappingMatrix = steeringMatrix;
    % -------------------------------------------------------------------
        
    % Generate a packet with random PSDU
    psduLength = getPSDULength(cfgHE); % PSDU length in bytes
    txPSDU = cell(numUsers,1);
    for userIdxIter = 1:numUsers
        txPSDU{userIdxIter} = randi([0 1],psduLength(userIdxIter)*8,1,'int8');
    end
    
    tx = wlanWaveformGenerator(txPSDU,cfgHE);
    % Add trailing zeros to allow for channel delay
    txPad = [tx; zeros(50,cfgHE.NumTransmitAntennas)];
    % Passing txPad into tgaxChannel
    [rx,pathGains] = tgaxChannel{userIdx}(txPad);

    % Get perfect timing offset and channel matrix for HE-LTF field
    heltfPathGains = pathGains(ind.HELTF(1):ind.HELTF(2),:,:,:,:);
    pktOffset = channelDelay(heltfPathGains,pathFilters);
    chan = helperPerfectChannelEstimate(heltfPathGains,pathFilters,ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices,pktOffset);

    % Pass the waveform through AWGN channel
    rx = awgnChannel(rx);

    % Demodulate data symbols
    rxData = rx(pktOffset+(ind.HEData(1):ind.HEData(2)),:);
    demodSym = wlanHEDemodulate(rxData,'HE-Data',cfgHE,ruIdx);

    % Extract data subcarriers from demodulated symbols and channel
    % estimate
    demodDataSym = demodSym(ofdmInfo.DataIndices,:,:);

    % Get channel estimate from channel matrix (include spatial mapping
    % and cyclic shift)
    chanEst = heChannelToChannelEstimate(chan,cfgHE,ruIdx);
    chanEstAv = permute(mean(chanEst,2),[1 3 4 2]); % Average over symbols
    chanEstData = chanEstAv(ofdmInfo.DataIndices,:,:);

    % Calculate single stream pilot estimates per symbol and noise
    % estimate
    chanEstSSPilots = permute(sum(chanEst(ofdmInfo.PilotIndices,:,:,:),3),[1 2 4 5 3]);
    demodPilotSym = demodSym(ofdmInfo.PilotIndices,:,:);
    nVarEst = heNoiseEstimate(demodPilotSym,chanEstSSPilots,cfgHE,ruIdx);

    % Equalization and STBC combining
    [eqDataSym,csi] = heEqualizeCombine(demodDataSym,chanEstData,nVarEst,cfgHE,userIdx);
    rxPSDU = wlanHEDataBitRecover(eqDataSym,nVarEst,csi,cfgHE,userIdx);

    % Determine if any bits are in error, i.e. a packet error
    packetError = ~isequal(txPSDU{userIdx},rxPSDU);
    perStore(numPkt) = packetError;
    numPacketErrors = numPacketErrors+packetError;

    numPkt = numPkt+1;
end

% Remove last increment
numPkt = numPkt-1;

% Calculate packet error rate (PER) at SNR point
packetErrorRate = numPacketErrors/numPkt;

% Return results
out = struct;
out.packetErrorRate = packetErrorRate;
out.perStore = perStore;
out.numPkt = numPkt;

disp([char(simParams.DelayProfile) ' '...
      num2str(simParams.NumTransmitAntennas) '-by-' ...
      num2str(simParams.NumReceiveAntennas) ','...
      ' MCS ' num2str(simParams.MCS) ','...
      ' SNR ' num2str(simParams.SNR) ...
      ' completed after ' num2str(out.numPkt) ' packets,'...
      ' PER:' num2str(out.packetErrorRate)]);
  
end