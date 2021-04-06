% Full PHY simulation of interferer-desired RX pair
% under 20MHz OFDM/OFDMA MIMO/MU-MIMO system
function out = box0Simulation(simParams)
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
tgaxChannel2 = cell(1,numUsers);
% Generate per-user channels
for userIdxIter = 1:numUsers
    tgaxChannel2{userIdxIter} = clone(simParams.Channel);
    tgaxChannel2{userIdxIter}.UserIndex = userIdxIter; % Set unique user index
end
tgaxChannelInterference = clone(simParams.Channel);

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
chanInfo = info(tgaxChannelInterference);
pathFilters = chanInfo.ChannelFilterCoefficients; % [NP numChannelTaps]

% Create object to deal with abstraction
Abstraction = tgaxLinkPerformanceModel;

% Loop to simulate multiple packets
perStore = nan(maxNumPackets,1);
perAbsStore = nan(maxNumPackets,1);
perAbsRawStore = nan(maxNumPackets,1);
snreffStore = nan(maxNumPackets,1);
sinrStore = nan(ruInfo(cfgHE).RUSizes(ruIdx),cfgHE.User{userIdx}.NumSpaceTimeStreams,maxNumPackets); % Nsc-by-Nsts-by-maxNumPackets
numPacketErrors = 0;
numPacketErrorsAbs = 0;
numPkt = 1; % Index of packet transmitted
while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
    % Reset BSS 2 channels for different realization
    for userIdxIter = 1:numel(tgaxChannel2)
        reset(tgaxChannel2{userIdxIter}); 
    end
    % Reset channel between interferer-desired RX pair for different realization
    reset(tgaxChannelInterference); 
    
    % The below dashed part caluclates ZF precoding matrix  
    % Comment the below dashed part to run simulation with default Fourier precoding
    % -------------------------------------------------------------------
    % Generate NDP packet - with an empty PSDU as no data
    txNDP = wlanWaveformGenerator([],cfgNDP);
    % For each user STA, pass the NDP packet through the channel and calculate
    % the feedback channel state matrix by SVD
    staFeedback = cell(1,numUsers);
    for userIdxIter = 1:numel(tgaxChannel2)
        % Received waveform at user STA with 50 sample padding. No noise.
        rxNDP = tgaxChannel2{userIdxIter}([txNDP; zeros(50,size(txNDP,2))]);
        
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
    [rx,pathGains] = tgaxChannelInterference(txPad);

    % Get perfect timing offset and channel matrix for HE-LTF field
    heltfPathGains = pathGains(ind.HELTF(1):ind.HELTF(2),:,:,:,:);
    pktOffset = channelDelay(heltfPathGains,pathFilters);
    chan = helperPerfectChannelEstimate(heltfPathGains,pathFilters,ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices,pktOffset);
    
    % Calculate SINR using abstraction
    % As multiple symbols returned average over symbols and permute
    % for calculations
    % Get precoding matrix for abstraction
    Wtx = getPrecodingMatrix(cfgHE,ruIdx); % Include cyclic shift applied per STS
    WtxUser1 = Wtx(:,1:Nsts,:); % Hard code
    WtxUser1 = WtxUser1/sqrt(Nsts);
    Htxrx = permute(mean(chan,2),[1 3 4 2]); % Nst-by-Nt-by-Nr
    numUsersCurrentRU = length(cfgHE.RU{ruIdx}.UserNumbers);
    Ptxrx = 1/numUsersCurrentRU; % Assume transmit power is 0dBW and uniformly splitted to each user
    sinr = calculateSINR(Htxrx,Ptxrx,WtxUser1,N0);
    sinrStore(:,:,numPkt) = sinr;
    
    % Link performance model - estimate PER using abstraction
    [perAbs,effSINR] = estimateLinkPerformance(Abstraction,sinr,cfgHE,userIdx);

    % Flip a coin for the abstracted PHY
    packetErrorAbs = rand(1)<=perAbs;
    numPacketErrorsAbs = numPacketErrorsAbs+packetErrorAbs;

    % Store outputs for analysis
    perAbsRawStore(numPkt) = perAbs;
    perAbsStore(numPkt) = packetErrorAbs;
    snreffStore(numPkt) = effSINR;

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

% Calculate packet error rate (PER) at SNR point
packetErrorRate = numPacketErrors/numPkt;
packetErrorRateAbs = numPacketErrorsAbs/numPkt;

% Return results
out = struct;
out.packetErrorRate = packetErrorRate;
out.perStore = perStore;
out.numPkt = numPkt;
out.sinrStore = sinrStore;
out.packetErrorRateAbs = packetErrorRateAbs;
out.perAbsRawStore = perAbsRawStore;
out.perAbsStore = perAbsStore;

disp([char(simParams.DelayProfile) ' '...
      num2str(simParams.NumTransmitAntennas) '-by-' ...
      num2str(simParams.NumReceiveAntennas) ','...
      ' MCS ' num2str(simParams.MCS) ','...
      ' SNR ' num2str(simParams.SNR) ...
      ' completed after ' num2str(out.numPkt) ' packets,'...
      ' PER:' num2str(out.packetErrorRate)]);
  
end