function steeringMat = heDirectUserBeamformingFeedback(tgax,cfgNDP,varargin)
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
% Indices of active subcarriers in the NDP
ndpOFDMInfo = wlanHEOFDMInfo('HE-Data',cfgNDP);

pathGains = tgax();
chanInfo = info(tgax);
pathFilters = chanInfo.ChannelFilterCoefficients; % [NP numChannelTaps]

chan = helperPerfectChannelEstimate(pathGains,pathFilters,ndpOFDMInfo.FFTLength,ndpOFDMInfo.CPLength,ndpOFDMInfo.ActiveFFTIndices);
chanEst = heChannelToChannelEstimate(chan,cfgNDP);
chanEstAv = permute(mean(chanEst,2),[1 3 4 2]); % Average over symbols

% Get cyclic shift and inverse
numSTS = cfgNDP.NumSpaceTimeStreams;
csh = -wlan.internal.getCyclicShiftVal('VHT',numSTS,wlan.internal.cbwStr2Num(chanBW));
Nfft = (fs/20e6)*256;
k = ndpOFDMInfo.ActiveFrequencyIndices;
chanEstMinusCSD = permute(wlan.internal.cyclicShift(permute(chanEstAv,[1 3 2]),csh,Nfft,k),[1 3 2]);
chanEstPerm = permute(chanEstMinusCSD,[3 2 1]); % Nr-by-Nsts-by-Nst

% Compute the feedback matrix using singular value decomposition
% for the streams allocated to the user
Nst = size(chanEstPerm,3);
numRx = size(chanEstPerm,1);
V = complex(zeros(Nst,numSTS,numRx)); % Nst-by-Nsts-by-Nr
for ist = 1:Nst
    [~,~,V(ist,:,:)] = svd(chanEstPerm(:,:,ist),'econ');
end
steeringMat = permute(V,[1 3 2]); % Nst-by-Nr-by-Nsts

end
