function steeringMatBF = heSUCalculateSteeringMatrix(steeringMatFB,cfg,cfgNDP)
%  heSUCalculateSteeringMatrix Calculate beamforming steering matrix
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   STEERINGMATBF = heSUCalculateSteeringMatrix(STEERINGMATFB,CFG,CFGNDP) 
%   returns the steering matrix recommended to beamform an iser in a transmit
%   beamforming 
%
%   STEERINGMATFB is a cell array containing the steering matrices fed-back
%   by each user in the RU to beamform.
%
%   CFG is the configuration of the HE-SU transmission
%   
%   CFGNDP is the configuration of the HE-NDP used to gather feedback and
%   is a format configuration object of type
%   <a href="matlab:help('wlanHESUConfig')">wlanHESUConfig</a>.
%

%   Copyright 2018 The MathWorks, Inc.

allocInfo = ruInfo(cfg);

% Indices of active subcarriers within the RU
OFDMInfo = wlanHEOFDMInfo('HE-Data',cfg);
ofdmInd = OFDMInfo.ActiveFrequencyIndices;

% Indices of active subcarriers in the NDP
ndpOFDMInfo = wlanHEOFDMInfo('HE-Data',cfgNDP);
trainingInd = ndpOFDMInfo.ActiveFrequencyIndices;

% Get the indices which overlap - use to extract from NDP
[~,scUseInd] = intersect(trainingInd,ofdmInd);

% Only take the RU subcarriers and space-time streams of
% interest for the current RU and user
numSTS = cfg.NumSpaceTimeStreams;
numRx = size(steeringMatFB,2);
if numSTS>numRx
    error('The number of space-time streams (%d) exceeds the number of receive antennas (%d) for user %d',numSTS,numRx);
end
steeringMatBF = steeringMatFB(scUseInd,1:numSTS,:);

end