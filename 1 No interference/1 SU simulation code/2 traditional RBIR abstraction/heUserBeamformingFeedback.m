function steeringMat = heUserBeamformingFeedback(rx,cfgNDP,varargin)
%heUserBeamformingFeedback HE user beamforming feedback
%
%   STEERINGMAT = heUserBeamformingFeedback(RX,CFGNDP) returns the steering
%   matrix recommended to beamform towards the user of interest. The
%   steering matrix is calculated using SVD.
%
%   STEERINGMAT = heUserBeamformingFeedback(RX,CFGNDP,CFOCOMP) performs
%   carrier frequency offset estimation and compensation on the received
%   waveform before calculating the steering matrix if CFOCOMP is set to
%   TRUE.
%
%   STEERINGMAT is a Nst-by-Nr-by-Nsts array containing the recommended
%   full band beamforming steering matrix. Nst is the number of occupied
%   subcarriers, Nr is the number of receive antennas, and Nsts is the
%   number of space-time-streams.
%
%   RX is the received NDP at the station.
%
%   CFGNDP is the format configuration object of type <a href="matlab:help(wlanHESUConfig')">wlanHESUConfig</a>.

%   Copyright 2017-2018 The MathWorks, Inc.

narginchk(2,3);

if nargin == 3
    cfoCompensate = varargin{1};
    validateattributes(cfoCompensate,{'logical'},{'nonnan','finite'},mfilename,'',3);
else
    cfoCompensate = false;
end

chanBW = string(cfgNDP.ChannelBandwidth);
ind = wlanFieldIndices(cfgNDP);
fs = wlanSampleRate(cfgNDP);

% Packet detect and determine coarse packet offset
coarsePktOffset = wlanPacketDetect(rx,chanBW);
if isempty(coarsePktOffset) % If empty no L-STF detected; try 0
    % Synchronization failed, return empty steering matrix
    steeringMat = [];
    return;
end

if cfoCompensate
    % Extract L-STF and perform coarse frequency offset correction
    lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
    coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
    rx = helperFrequencyOffset(rx,fs,-coarseFreqOff);
end

% Extract the non-HT fields and determine fine packet offset
nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);

% Determine final packet offset
pktOffset = coarsePktOffset+finePktOffset;

% If packet detected outwith the range of expected delays from
% the channel modeling; packet error
if pktOffset>50
    % Synchronization failed, return empty steering matrix
    steeringMat = [];
    return;
end

if cfoCompensate
    % Extract L-LTF and perform fine frequency offset correction
    rxLLTF = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
    fineFreqOff = wlanFineCFOEstimate(rxLLTF,chanBW);
    rx = helperFrequencyOffset(rx,fs,-fineFreqOff);
end

% Demodulate HE-LTF with info for all RUs
rxHELTF = rx(pktOffset+(ind.HELTF(1):ind.HELTF(2)),:);
demodHELTF = wlanHEDemodulate(rxHELTF,'HE-LTF',cfgNDP);

% Extract the demodulated HE-LTF RU
demodHELTFRU = demodHELTF;

% Channel estimate with info from current RU
chanEstUser = heLTFChannelEstimate(demodHELTFRU,cfgNDP);

% Get cyclic shift and inverse
numSTS = cfgNDP.NumSpaceTimeStreams;
csh = -wlan.internal.getCyclicShiftVal('VHT',numSTS,wlan.internal.cbwStr2Num(chanBW));
Nfft = (fs/20e6)*256;
% Indices of active subcarriers in the NDP
ndpOFDMInfo = wlanHEOFDMInfo('HE-Data',cfgNDP);
k = ndpOFDMInfo.ActiveFrequencyIndices;
chanEstMinusCSD = permute(wlan.internal.cyclicShift(permute(chanEstUser,[1 3 2]),csh,Nfft,k),[1 3 2]);
chanEstPerm = permute(chanEstMinusCSD,[3 2 1]); % Nr-by-Nsts-by-Nst

% Compute the feedback matrix using singular value decomposition
% for the streams allocated to the user
Nst = size(chanEstPerm,3);
numRx = size(demodHELTF,3);
V = complex(zeros(Nst,numSTS,numRx)); % Nst-by-Nsts-by-Nr
for ist = 1:Nst
    [~,~,V(ist,:,:)] = svd(chanEstPerm(:,:,ist),'econ');
end
steeringMat = permute(V,[1 3 2]); % Nst-by-Nr-by-Nsts

end
