function [chanEstRU,varargout] = heLTFChannelEstimate(demodHELTFRU,cfg,varargin)
%heLTFChannelEstimate Channel estimation using HE-LTF
%   CHANESTRU = heLTFChannelEstimate(RXSYM,CFGHE) returns the estimated
%   channel between all space-time streams and receive antennas using
%   HE-LTF of an HE single user, extended range single user, multi-user or
%   trigger-based (HE TB) packet. The channel estimate includes the
%   effect of the applied spatial mapping matrix and cyclic shifts at the
%   transmitter. If HE-LTF compression is used, linear interpolation is
%   performed to create a channel estimate for all subcarriers.
%
%   CHANESTRU is an array characterizing the estimated channel for the data
%   and pilot subcarriers. EST is a complex Nst-by-Nsts-by-Nr array
%   characterizing the estimated channel for the data and pilot
%   subcarriers, where Nst is the number of occupied subcarriers, Nsts is
%   the total number of space-time streams, and Nr is the number of receive
%   antennas. If CFGHE is a MU configuration, then the channel estimate for
%   all RUs is returned.
%
%   RXSYM is a complex Nst-by-Nsym-by-Nr array containing demodulated
%   concatenated HE-LTF. Nsym is the number of demodulated HE-LTF symbols.
%
%   CFGHE is a format configuration object of type <a href="matlab:help('wlanHESUConfig')">wlanHESUConfig</a>
%   <a href="matlab:help('wlanHEMUConfig')">wlanHEMUConfig</a>, <a href="matlab:help('wlanHETBConfig')">wlanHETBConfig</a>, <a href="matlab:help('heTBSystemConfig')">heTBSystemConfig</a> or <a href="matlab:help('wlanHERecoveryConfig')">wlanHERecoveryConfig</a>.
%
%   CHANESTRU = heLTFChannelEstimate(RXSYM,CFGMU,RUOFINTEREST) returns the
%   channel estimate for the RU of interest index RUOFINTEREST for a
%   multi-user configuration. CFGMU is of type <a href="matlab:help('wlanHEMUConfig')">wlanHEMUConfig</a> or 
%   <a href="matlab:help('heTBSystemConfig')">heTBSystemConfig</a>. If not provided the default is 1.
%
%   [...,CHANESTSSPILOTS] = heLTFChannelEstimate(...) additionally returns
%   a Nsp-by-Nsym-by-Nr array characterizing the estimated channel for
%   pilot subcarrier locations for each symbol, assuming one space-time
%   stream at the transmitter.
%
%   Examples:
%   % Decode the HE data field for each user in and OFDMA and MU-MIMO
%   % transmission with a fading channel model. Estimate the channel for
%   % each user.
%
%      % Create packet configuration
%      allocationIndex = [192 193]; % Two 242 RUs, the second with 2 users
%      cfg = wlanHEMUConfig(allocationIndex);
%      cfg.NumTransmitAntennas = 2;
%      cfg.User{1}.NumSpaceTimeStreams = 2;
% 
%      % Generate MU waveform
%      txWaveform = wlanWaveformGenerator([1;0;0;1],cfg);
% 
%      % Channel and receiver per user
%      ind = wlanFieldIndices(cfg);
%      allocationInfo = ruInfo(cfg);
%      for ruIdx = 1:allocationInfo.NumRUs
%          for userIdx = 1:allocationInfo.NumUsersPerRU(ruIdx)
%              % Add channel and noise
%              snr = 20;
%              channel = wlanTGaxChannel;
%              channel.NumTransmitAntennas = 2;
%              channel.NumReceiveAntennas = 2;
%              channel.SampleRate = wlanSampleRate(cfg);
%              rxWaveform = awgn(channel([txWaveform; zeros(10,2)]),snr);
% 
%              % Synchronize
%              rxWaveform = rxWaveform(1+4:end,:);
% 
%              % Extract and OFDM demodulate the HE-LTF for the RU of
%              % interest
%              rxHETLF = rxWaveform(ind.HELTF(1): ind.HELTF(2),:);
%              demodHELTF = wlanHEDemodulate(rxHETLF,'HE-LTF',cfg,ruIdx);
% 
%              % Channel estimate for RU of interest
%              chanEst = heLTFChannelEstimate(demodHELTF,cfg,ruIdx);
% 
%              % Extract and OFDM demodulate the data field for the RU of
%              % interest
%              rxData = rxWaveform(ind.HEData(1):ind.HEData(2),:);
%              demodData = wlanHEDemodulate(rxData,'HE-Data',cfg,ruIdx);
% 
%              % Equalize data symbols - extract the space-time streams for
%              % the user of interest after equalization
%              nVar = 10^-(snr/10);
%              [eqSym,csi] = heEqualizeCombine(demodData,chanEst,nVar, ...
%                 cfg,userIdx);
% 
%              % Discard pilot carriers and decode
%              info = wlanHEOFDMInfo('HE-Data',cfg,ruIdx);
%              rxBits = wlanHEDataBitRecover(eqSym(info.DataIndices,:,:), ...
%                  nVar,csi(info.DataIndices,:),cfg,userIdx);
%          end
%      end
 
%   Copyright 2017-2019 The MathWorks, Inc.

%#codegen

if nargin>2
    ruOfInterest = varargin{1};
else
    ruOfInterest = 1;
end

% Validate the format configuration object is a valid type
validateattributes(cfg,{'wlanHESUConfig','wlanHEMUConfig','wlanHETBConfig','heTBSystemConfig','wlanHERecoveryConfig'},{'scalar'},mfilename,'format configuration object');

% Get allocation information
if isa(cfg,'wlanHERecoveryConfig')
    ruSizeRU = cfg.RUSize;
    ruIndexRU = cfg.RUIndex;
    pktFormat = cfg.PacketFormat;
    if strcmp(pktFormat,'HE-MU')
        numSTSRU = cfg.RUTotalSpaceTimeStreams;
    else % SU, EXT SU
        numSTSRU = cfg.NumSpaceTimeStreams;
    end
else
    allocInfo = ruInfo(cfg);
    coder.internal.errorIf(ruOfInterest>allocInfo.NumRUs,'wlan:he:InvalidRUOfInterest',ruOfInterest,allocInfo.NumRUs);
    ruSizeRU = allocInfo.RUSizes(ruOfInterest);
    ruIndexRU = allocInfo.RUIndices(ruOfInterest);
    numSTSRU = allocInfo.NumSpaceTimeStreamsPerRU(ruOfInterest);
    pktFormat = packetFormat(cfg);
end

% Validate symbol type
validateattributes(demodHELTFRU,{'single','double'},{'3d'},mfilename,'HE-LTF OFDM symbol(s)');
[numST,numLTF,numRx] = size(demodHELTFRU);
tac = wlan.internal.heRUToneAllocationConstants(ruSizeRU);
coder.internal.errorIf(numST~=tac.NST,'wlan:wlanChannelEstimate:IncorrectNumSC',tac.NST,numST);
ofdmInfo = wlanHEOFDMInfo('HE-LTF',cfg.ChannelBandwidth,cfg.GuardInterval,[ruSizeRU ruIndexRU]);
if numLTF==0
    chanEstRU = zeros(numST,numSTSRU,numRx);
    varargout{1} = zeros(numel(ofdmInfo.PilotIndices),numLTF,numRx); % For codegen
    return;
end
minNumLTF = wlan.internal.numVHTLTFSymbols(numSTSRU);
coder.internal.errorIf(numLTF<minNumLTF,'wlan:he:InvalidNumLTF',numLTF,minNumLTF);

% Get the HE-LTF sequence
cbw = wlan.internal.cbwStr2Num(cfg.ChannelBandwidth);
[HELTF,kHELTFSeq] = wlan.internal.heLTFSequence(cbw,cfg.HELTFType);

% Extract the RU of interest from the full-bandwidth HELTF
kRU = ofdmInfo.ActiveFrequencyIndices;
[~,ruIdx] = intersect(kHELTFSeq,kRU);
HELTFRU = HELTF(ruIdx);

switch cfg.HELTFType
    % IEEE P802.11ax/D4.1, Equation 27-52
    case 1
        N_HE_LTF_Mode = 4; % undefined
    case 2
        N_HE_LTF_Mode = 2;
    otherwise % 4
        N_HE_LTF_Mode = 1;
end

isaTBConfig = isa(cfg,'heTBSystemConfig') || isa(cfg,'wlanHETBConfig');
if numSTSRU==1
    % Single STS

    % When more than one LTF we can average over the LTFs for data and
    % pilots to improve the estimate. As there is only one space-time
    % stream, the pilots and data essentially both use the P matrix which
    % does not change per space-time stream (only per symbol), therefore
    % this "MIMO" estimate performs the averaging of the number of symbols.
    chanEstRU = wlan.internal.mimoChannelEstimate(demodHELTFRU,HELTFRU,numSTSRU);
    
    % Remove orthogonal sequence across subcarriers (if used)
    if isaTBConfig && cfg.SingleStreamPilots==false
        chanEstRU = removeOrthogonalSequence(chanEstRU,numSTSRU,kRU,N_HE_LTF_Mode);
    end

    % Interpolate if HE-LTF compression used
    if N_HE_LTF_Mode>1
        chanEstRU = chanEstInterp(chanEstRU,cbw,N_HE_LTF_Mode,ruSizeRU,ruIndexRU);
    end
else 
    % MIMO channel estimation as per Perahia, Eldad, and Robert Stacey.
    % Next Generation Wireless LANs: 802.11 n and 802.11 ac. Cambridge
    % University Press, 2013, page 100, Equation 4.39.
    % Remove orthogonal sequence across subcarriers (if used)
    if isaTBConfig && cfg.SingleStreamPilots==false
        % Only perform channel estimate for non-pilot subcarriers as pilots
        % are single stream
        kMIMO = kRU; % All subcarriers MIMO estimates 
        mimoInd = 1:numST;
        chanEstRUMIMO = wlan.internal.mimoChannelEstimate(demodHELTFRU,HELTFRU,numSTSRU);
        chanEstRUMIMO = removeOrthogonalSequence(chanEstRUMIMO,numSTSRU,kRU,N_HE_LTF_Mode);
    else
        % Only perform channel estimate for non-pilot subcarriers as pilots
        % are single stream
        mimoInd = ofdmInfo.DataIndices;
        kMIMO = kRU(mimoInd); % Only data subcarriers MIMO estimates
        chanEstRUMIMO = wlan.internal.mimoChannelEstimate(demodHELTFRU(ofdmInfo.DataIndices,:,:),HELTFRU(mimoInd),numSTSRU);
    end
    
    % Undo cyclic shift for each STS before averaging and interpolation    
    nfft = (cbw/20)*256;
    numSTSTotal = size(chanEstRUMIMO,2);
    csh = wlan.internal.getCyclicShiftVal('VHT',numSTSTotal,cbw);
    chanEstRUMIMO = wlan.internal.cyclicShiftChannelEstimate(chanEstRUMIMO,-csh,nfft,kMIMO);

    % Interpolate over pilot locations, and any compressed subcarriers
    chanEstRU = chanEstInterp(chanEstRUMIMO,cbw,N_HE_LTF_Mode,ruSizeRU,ruIndexRU,mimoInd);

    % Re-apply cyclic shift after interpolation
    chanEstRU = wlan.internal.cyclicShiftChannelEstimate(chanEstRU,csh,nfft,kRU);
end

% If extended range SU, then the HE-LTF are boosted by sqrt(2). If we
% don't remove this at demodulation then we must de-scale the channel
% estimate as the data field is not scaled.
if strcmp(pktFormat,'HE-EXT-SU')
    eta = 1/sqrt(2);
else
    eta = 1;
end
chanEstRU = chanEstRU*eta; % Scale for HE-EXT-SU

% Channel estimate for pilots
if nargout>1
    if isaTBConfig && cfg.SingleStreamPilots==false
        % Create single stream from MIMO pilot estimates by summing across
        % space-time streams (2nd dimension)
        varargout{1} = sum(chanEstRU(ofdmInfo.PilotIndices,:,:),2);
    else
        % Channel estimate for single-stream pilots
        Pheltf = wlan.internal.mappingMatrix(numLTF);
        R = Pheltf(1,1:numLTF); % R matrix changes pilot polarity per symbol
        % Estimate the channel at pilot subcarriers accounting for polarity
        chanEstSSPilots = bsxfun(@rdivide,demodHELTFRU(ofdmInfo.PilotIndices,:,:),bsxfun(@times,HELTFRU(ofdmInfo.PilotIndices),R));
        varargout{1} = chanEstSSPilots*eta; % Scale for HE_EXT_SU
    end
end

end

function chanEstRUInterp = chanEstInterp(chanEstRU,cbw,N_HE_LTF_Mode,ruSize,ruIndex,varargin)
    % Interpolate over pilot locations and compressed subcarriers

    Nfft = 256*cbw/20;
    
    % Get the subcarrier indices within the FFT for the channel estimate
    % input
    kAct = wlan.internal.heRUSubcarrierIndices(cbw,ruSize,ruIndex)+Nfft/2+1;
    % If the channelEstRU is not the entire RU, then we need to make sure
    % we know the subcarrier indices, so use the ruInd input. For example
    % this allows us to pass in only the data subcarriers.
    if nargin>5
        ruInd = varargin{1};
        kChanEstInputs = kAct(ruInd);
    else
        % Assume chanEstRU is the whole RU
        kChanEstInputs = kAct;
    end
    
    % Get the indices within the FFT which contain actual estimates
    % (excluding the guard bands). This is how the pattern is structured
    kAll = 1:N_HE_LTF_Mode:Nfft; 
    
    % Find the subcarrier indices within the FFT which contain actual data
    % within the channel estimate input (kToInterp) and the indices of
    % these within the chanEstDataRU input array (toInterpInd)
    [kToInterp,toInterpInd] = intersect(kChanEstInputs,kAll);

    % Interpolate and extrapolate over all RU subcarrier indices to
    % interpolate over compressed region and pilots
    magPart = interp1(kToInterp.',abs(chanEstRU(toInterpInd,:,:)),kAct,'linear','extrap');
    phasePart = interp1(kToInterp.',unwrap(angle(chanEstRU(toInterpInd,:,:))),kAct,'linear','extrap');
    [realPart,imagPart] = pol2cart(phasePart,magPart);
    chanEstRUInterp = complex(realPart,imagPart);

end

function chanEstRUData = removeOrthogonalSequence(chanEstRUData,numSTSRU,k,N_HE_LTF_Mode)
    % Remove the orthogonal sequence across subcarriers
    M = 0; % Assume space-time streams of all users in estimate
    m = 1:numSTSRU;
    Pheltf = wlan.internal.mappingMatrix(8);
    seq = Pheltf(M+m,mod(ceil(k/N_HE_LTF_Mode)-1,8)+1).'; % Nsts-by-Nst
    chanEstRUData = chanEstRUData./seq;
end
