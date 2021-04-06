% Full PHY simulation with interference
% under two 20MHz OFDM/OFDMA MIMO/MU-MIMO systems with same setup
function out = box0Simulation(simParams)
% box0Simulation Example helper function

% Copyright 2019 The MathWorks, Inc.

% Extract configuration
userIdx = simParams.userIdx;
ruIdx = simParams.ruIdx;
cfgHE = simParams.Config;
cfgHEInt = simParams.Config;
intPathlossdB = simParams.intPathlossdB;
substreamidx = simParams.RandomSubstream;
maxNumPackets = simParams.MaxNumPackets;
maxNumErrors = simParams.MaxNumErrors;
numUsers = size(cfgHE.User,2);
snr = simParams.SNR;
% Generate per-user channels for BSS 1 (victim BSS)
tgaxChannel = cell(1,numUsers);
for userIdxIter = 1:numUsers
    tgaxChannel{userIdxIter} = clone(simParams.Channel);
    tgaxChannel{userIdxIter}.UserIndex = userIdxIter; % Set unique user index
end
% Generate per-user channels for BSS 2 (interfering BSS)
tgaxChannel2 = cell(1,numUsers);
for userIdxIter = 1:numUsers
    tgaxChannel2{userIdxIter} = clone(simParams.Channel);
    tgaxChannel2{userIdxIter}.UserIndex = userIdxIter; % Set unique user index
end
% Set up interference channel from BSS 2 to BSS 1
tgaxChannelInterference = clone(simParams.Channel);

% Create an NDP packet for BSS 1
cfgNDP = wlanHESUConfig('APEPLength',0,'GuardInterval',0.8); % No data in an NDP
cfgNDP.ChannelBandwidth = cfgHE.ChannelBandwidth;
cfgNDP.NumTransmitAntennas = cfgHE.NumTransmitAntennas;
cfgNDP.NumSpaceTimeStreams = cfgHE.NumTransmitAntennas;

% Create an NDP packet for BSS 2 
cfgNDPInt = wlanHESUConfig('APEPLength',0,'GuardInterval',0.8); % No data in an NDP
cfgNDPInt.ChannelBandwidth = cfgHEInt.ChannelBandwidth;
cfgNDPInt.NumTransmitAntennas = cfgHEInt.NumTransmitAntennas;
cfgNDPInt.NumSpaceTimeStreams = cfgHEInt.NumTransmitAntennas;

% Set random substream index per iteration to ensure that each
% iteration uses a repeatable set of random numbers
stream = RandStream('combRecursive','Seed',99);
stream.Substream = substreamidx;
RandStream.setGlobalStream(stream);

% Indices to extract fields from the PPDU
ind = wlanFieldIndices(cfgHE);
indInt = wlanFieldIndices(cfgHEInt);

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
pathFilters = chanInfo.ChannelFilterCoefficients; 
interferenceChanInfo = info(tgaxChannelInterference);
interferencePathFilter = interferenceChanInfo.ChannelFilterCoefficients;

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
    % Reset channel for different realization
    for userIdxIter = 1:numel(tgaxChannel)
        reset(tgaxChannel{userIdxIter}); % Reset BSS 1 channels 
        reset(tgaxChannel2{userIdxIter}); % Reset BSS 2 channels 
    end
    % Reset interference channel for different realization
    reset(tgaxChannelInterference);
    
    % The below dashed part caluclates ZF precoding matrix for BSS 1
    % Comment the below dashed part to run simulation with default Fourier precoding
    % -------------------------------------------------------------------
    % Generate NDP packet in BSS 1 - with an empty PSDU as no data
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

    % The below dashed part caluclates ZF precoding matrix for BSS 2
    % Comment the below dashed part to run simulation with default Fourier precoding
    % -------------------------------------------------------------------
    % Generate NDP packet - with an empty PSDU as no data
    txNDPInt = wlanWaveformGenerator([],cfgNDPInt);
    % For each user STA, pass the NDP packet through the channel and calculate
    % the feedback channel state matrix by SVD
    staFeedbackInt = cell(1,numUsers);
    for userIdxIter = 1:numel(tgaxChannel2)
        % Received waveform at user STA with 50 sample padding. No noise.
        rxNDPInt = tgaxChannel2{userIdxIter}([txNDPInt; zeros(50,size(txNDPInt,2))]);
        
        % Get the full-band beamforming feedback for a user
        staFeedbackInt{userIdxIter} = heUserBeamformingFeedback(rxNDPInt,cfgNDPInt);
    end
    % Calculate the steering matrix to apply to the RU given the feedback
    % A zero forcing solution is used to calculate the steering matrix
    steeringMatrixInt = heMUCalculateSteeringMatrix(staFeedbackInt,cfgHEInt,cfgNDPInt,ruIdx);
    
    % Apply the steering matrix to the RU of interest
    cfgHEInt.RU{ruIdx}.SpatialMapping = 'Custom';
    cfgHEInt.RU{ruIdx}.SpatialMappingMatrix = steeringMatrixInt;
    % -------------------------------------------------------------------
    
    % Setup of the desired signal in BSS 1
    psduLength = getPSDULength(cfgHE); % PSDU length in bytes
    txPSDU = cell(numUsers,1);
    for userIdxIter = 1:numUsers
        txPSDU{userIdxIter} = randi([0 1],psduLength(userIdxIter)*8,1,'int8');
    end
    % Transmitted desired signal
    tx = wlanWaveformGenerator(txPSDU,cfgHE);
    % Add trailing zeros to allow for channel delay
    txPad = [tx; zeros(50,cfgHE.NumTransmitAntennas)];
    % Pass txPad into tgaxChannel and get received desired signal
    [rx,pathGains] = tgaxChannel{userIdx}(txPad);
    % Get perfect timing offset and channel matrix for HE-LTF field
    heltfPathGains = pathGains(ind.HELTF(1):ind.HELTF(2),:,:,:,:);
    pktOffset = channelDelay(heltfPathGains,pathFilters);
    chan = helperPerfectChannelEstimate(heltfPathGains,pathFilters,ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices,pktOffset);
    
    % Setup of interference signal in BSS 2
    psduLengthInt = getPSDULength(cfgHEInt); % PSDU length in bytes
    interferencePSDU = cell(numUsers,1);
    for userIdxIter = 1:numUsers
        interferencePSDU{userIdxIter} = randi([0 1],psduLengthInt(userIdxIter)*8,1,'int8');
    end
    % Transmitted interference signal
    txInt = wlanWaveformGenerator(interferencePSDU,cfgHEInt);
    % Add trailing zeros to allow for channel delay
    txIntPad = [txInt; zeros(50,cfgHEInt.NumTransmitAntennas)];
    % Pass txIntPad into tgaxChannelInterference and get received interference signal
    [rxInt,pathGainsInt] = tgaxChannelInterference(txIntPad);
    % Get perfect timing offset and channel matrix for HE-LTF field
    heltfPathGainsInt = pathGainsInt(indInt.HELTF(1):indInt.HELTF(2),:,:,:,:);
    pktOffsetInt = channelDelay(heltfPathGainsInt,interferencePathFilter);
    chanInt = helperPerfectChannelEstimate(heltfPathGainsInt,interferencePathFilter,ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices,pktOffsetInt);

    % Calculate SINR using abstraction
    % User 1: desired user
    % User 2: interfering source
    % Get precoding matrix for abstraction
    Wtx = getPrecodingMatrix(cfgHE,ruIdx); % Include cyclic shift applied per STS
    WtxUser1 = Wtx(:,1:Nsts,:);
    WtxUser1 = WtxUser1/sqrt(Nsts);
    WtxInt = getPrecodingMatrix(cfgHEInt,ruIdx); % Include cyclic shift applied per STS
    WtxUser2 = WtxInt(:,1:Nsts,:);
    WtxUser2 = WtxUser2/sqrt(Nsts); 
    % Get channel matrix for abstraction
    Htxrx = permute(mean(chan,2),[1 3 4 2]); % Nst-by-Nt-by-Nr
    Hint = permute(mean(chanInt,2),[1 3 4 2]); % Nst-by-Nt-by-Nr
    % Get transmit power for abstraction
    numUsersCurrentRU = length(cfgHE.RU{ruIdx}.UserNumbers);
    Ptxrx = 1/numUsersCurrentRU; % Assume transmit power is 0dBW and uniformly splitted to each user
    intPathloss =  1/10^(intPathlossdB/10); % Interference path loss in linear scale
    Pint = Ptxrx * intPathloss; % Applying interference path loss
    % Get post-MIMO processing SINR
    sinr = calculateSINR(Htxrx,Ptxrx,WtxUser1,N0,{Hint},Pint,{WtxUser2});
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

    % Assume interference packets and desired packets are aligned in time
    % before equalization
    rxData = rx(pktOffset+(ind.HEData(1):ind.HEData(2)),:);
    rxDataInt = rxInt(pktOffsetInt+(indInt.HEData(1):indInt.HEData(2)),:);
    
    % Pre-equalization combined desired signal and interference signal
    rxDataAdd = rxData + rxDataInt * sqrt(intPathloss);

    % Extract data subcarriers from demodulated symbols and channel
    % estimate
    demodSym = wlanHEDemodulate(rxDataAdd,'HE-Data',cfgHE,ruIdx);
    demodDataSym = demodSym(ofdmInfo.DataIndices,:,:);

    % Get channel estimate from channel matrix (include spatial mapping
    % and cyclic shift)
    chanEst = heChannelToChannelEstimate(chan,cfgHE,ruIdx); % Use desired channel estimate
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