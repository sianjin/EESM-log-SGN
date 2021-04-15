% Interference-free traditional EESM PHY layer abstraction
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
beta = simParams.beta;
% Generate per-user channels, adapted from HEDownlinkMUExample
tgaxChannel = cell(1,numUsers);
% Generate per-user channels
for userIdxIter = 1:numUsers
    tgaxChannel{userIdxIter} = clone(simParams.Channel);
    tgaxChannel{userIdxIter}.ChannelFiltering = false;
    tgaxChannel{userIdxIter}.NormalizePathGains = true;
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

% Create object to deal with abstraction
Abstraction = tgaxEESMLinkPerformanceModel;

% Loop to simulate multiple packets
numPacketErrors = 0;
numPacketErrorsAbs = 0;
numPkt = 1; % Index of packet transmitted
% Storing effective SNR vector for estimating effective SNR distribution
effSnrVec = [];
while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
    % Reset channel for different realization
    for userIdxIter = 1:numel(tgaxChannel)
        reset(tgaxChannel{userIdxIter}); 
    end
    
    % The below dashed part caluclates ZF precoding matrix  
    % Comment the below dashed part to run simulation with default Fourier precoding
    % -------------------------------------------------------------------
    % For each user STA, calculate the feedback channel state matrix by SVD
    staFeedback = cell(1,numUsers);
    for userIdxIter = 1:numel(tgaxChannel)     
        % Directly get the full-band beamforming feedback for a user
        staFeedback{userIdxIter} = heDirectUserBeamformingFeedback(tgaxChannel{userIdxIter},cfgNDP);
    end
    % Calculate the steering matrix to apply to the RU given the feedback
    % A zero forcing solution is used to calculate the steering matrix
    steeringMatrix = heMUCalculateSteeringMatrix(staFeedback,cfgHE,cfgNDP,ruIdx);
    
    % Apply the steering matrix to the RU of interest
    cfgHE.RU{ruIdx}.SpatialMapping = 'Custom';
    cfgHE.RU{ruIdx}.SpatialMappingMatrix = steeringMatrix;
    % -------------------------------------------------------------------
    
    pathGains = tgaxChannel{userIdx}();    % Get path gains
    chan = helperPerfectChannelEstimate(pathGains,pathFilters,ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices);
    
    % Calculate SINR using abstraction
    % As multiple symbols returned average over symbols and permute
    % for calculations
    % Get precoding matrix for abstraction
    Wtx = getPrecodingMatrix(cfgHE,ruIdx); % Include cyclic shift applied per STS
    WtxUser1 = Wtx(:,1:Nsts,:); % Hard code for the 1st user of current RU
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