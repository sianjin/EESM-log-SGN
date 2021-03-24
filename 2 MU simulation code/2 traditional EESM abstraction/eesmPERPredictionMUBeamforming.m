function out = eesmPERPredictionMUBeamforming(simParams)
% box0Simulation Example helper function

%   Copyright 2019 The MathWorks, Inc.

% Extract configuration
cfgHE = simParams.Config;
substreamidx = simParams.RandomSubstream;
maxNumPackets = simParams.MaxNumPackets;
maxNumErrors = simParams.MaxNumErrors;
numUsers = size(cfgHE.User,2);
snr = simParams.SNR;
beta = simParams.Beta;
% Generate per-user channels, adapted from HEDownlinkMUExample
tgaxChannel = cell(1,numUsers);
% Generate per-user channels
for userIdx = 1:numUsers
    tgaxChannel{userIdx} = clone(simParams.Channel);
    tgaxChannel{userIdx}.ChannelFiltering = false;
    tgaxChannel{userIdx}.NormalizePathGains = true;
    tgaxChannel{userIdx}.UserIndex = userIdx; % Set unique user index
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

% Per user processing, focus on the first user
userIdx = 1;
% Get occupied subcarrier indices and OFDM parameters
ofdmInfo = wlanHEOFDMInfo('HE-Data',cfgHE,userIdx);

% Create an instance of the AWGN channel per SNR point simulated
awgnChannel = comm.AWGNChannel;
awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
% Account for noise energy in nulls so the SNR is defined per
% active subcarrier
awgnChannel.SNR = snr-10*log10(ofdmInfo.FFTLength/sum(ruInfo(cfgHE).RUSizes));

% N0 = 10^(-snr/10); % original, I don't think it's correct
N0 = 10^(-awgnChannel.SNR/10); % This is correct

% Get path filers for last channel (same for all channels)
chanInfo = info(tgaxChannel{userIdx});
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
    % Pass through a fading indoor TGax channel
    for userIdx = 1:numel(tgaxChannel)
        reset(tgaxChannel{userIdx}); % Reset channel for different realization
    end
    
    % Per user processing, focus on the first user
    userIdx = 1;
    ruIdx = 1;
    
    pathGains = tgaxChannel{userIdx}();    % Get path gains
    chan = helperPerfectChannelEstimate(pathGains,pathFilters,ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices);
    
    % -------------------------------------------------------------------
    % For each user STA, pass the NDP packet through the channel and calculate
    % the feedback channel state matrix by SVD.
    staFeedback = cell(1,numUsers);
    for userIdx = 1:numel(tgaxChannel)
        staFeedback{userIdx} = heIdealUserBeamformingFeedback(tgaxChannel{userIdx},cfgNDP); % Get the full-band beamforming feedback for a user
    end
    steeringMatrix = heMUCalculateSteeringMatrix(staFeedback,cfgHE,cfgNDP,ruIdx);
    cfgHE.RU{ruIdx}.SpatialMapping = 'Custom';
    cfgHE.RU{ruIdx}.SpatialMappingMatrix = steeringMatrix;
        
%     % For each RU, calculate the steering matrix to apply
%     for ruIdx = 1:numel(cfgHE.RU)
%         % Calculate the steering matrix to apply to the RU given the feedback
%         steeringMatrix = heMUCalculateSteeringMatrix(staFeedback,cfgHE,cfgNDP,ruIdx);
%         
%         % Apply the steering matrix to each RU
%         cfgHE.RU{ruIdx}.SpatialMapping = 'Custom';
%         cfgHE.RU{ruIdx}.SpatialMappingMatrix = steeringMatrix;
%     end
    % -------------------------------------------------------------------
    
    % Per user processing, focus on the first user
    userIdx = 1;
    Nsts = cfgHE.User{userIdx}.NumSpaceTimeStreams;
    % Calculate SINR using abstraction
    % As multiple symbols returned average over symbols and permute
    % for calculations
    % Get precoding matrix for abstraction
    Wtx = getPrecodingMatrix(cfgHE,ruIdx); % Include cyclic shift applied per STS
    WtxUser1 = Wtx(:,1:Nsts,:);
    WtxUser1 = WtxUser1/sqrt(Nsts);
    Htxrx = permute(mean(chan,2),[1 3 4 2]); % Nst-by-Nt-by-Nr
    numUsersCurrentRU = length(cfgHE.RU{ruIdx}.UserNumbers);
    Ptxrx = 1/numUsersCurrentRU; % Assume transmit power is 0dBW and uniformly splitted to each user
    sinr = calculateSINR(Htxrx,Ptxrx,WtxUser1,N0);
    
    % Link performance model - estimate PER using abstraction
    snrEff = eesmEffectiveSINR(Abstraction,sinr,beta);
    effSnrVec = [effSnrVec, snrEff];
    perIns = estimateLinkPerformance(Abstraction,snrEff,cfgHE,userIdx);

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