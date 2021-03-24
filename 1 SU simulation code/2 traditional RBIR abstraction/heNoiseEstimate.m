function [nest,sigest] = heNoiseEstimate(x,chanEstSSPilots,cfg,varargin)
%heNoiseEstimate Estimate noise power using HE data field pilots
%
%   NEST = heNoiseEstimate(x,CHANESTSSPILOTS,CFGHE) estimates the mean
%   noise power in watts using the demodulated pilot symbols in the HE data
%   field and single-stream channel estimates at pilot subcarriers. The
%   noise estimate is averaged over the number of symbols and receive
%   antennas.
%
%   X is a complex Nsp-by-Nsym-by-Nr array containing demodulated pilot
%   subcarrier in HE data field. Nsym is the number of demodulated HE-Data
%   symbols.
%
%   CHANESTSSPILOTS is a complex Nsp-by-Nltf-by-Nr array containing the
%   channel gains at pilot subcarrier locations for each symbol, assuming
%   one space-time stream at the transmitter. Nltf is the number of HE-LTF
%   symbols.
%
%   CFGHE is a format configuration object of type <a href="matlab:help('wlanHESUConfig')">wlanHESUConfig</a>,
%   <a href="matlab:help('wlanHETBConfig')">wlanHETBConfig</a> or <a href="matlab:help('wlanHERecoveryConfig')">wlanHERecoveryConfig</a>.
%
%   NEST = heNoiseEstimate(X,CHANESTSSPILOTS,CFGMU,RUIDX) performs noise
%   power estimation for the multi user HE format input X.
%
%   CFGMU is the format configuration object of type <a href="matlab:help('wlanHEMUConfig')">wlanHEMUConfig</a> or
%   <a href="matlab:help('heTBSystemConfig')">heTBSystemConfig</a>.
%
%   RUIDX is the RU (resource unit) index.
%
%   [NEST,SIGEST] = heNoiseEstimate(...) additionally returns an estimate
%   of the signal power.

%   Copyright 2018-2019 The MathWorks, Inc.

%#codegen

validateattributes(cfg,{'wlanHESUConfig','wlanHEMUConfig','wlanHETBConfig','wlanHERecoveryConfig','heTBSystemConfig'},{'scalar'},mfilename,'format configuration object');

numOFDMSym = size(x,2);
n = (0:numOFDMSym-1);

ruIdx = 1;
if isa(cfg,'wlanHEMUConfig')
    narginchk(4,4)
    ruIdx = varargin{1};
    sigbInfo = wlan.internal.heSIGBCodingInfo(cfg);
    numHESIGB = sigbInfo.NumSymbols;
    pktFormat = packetFormat(cfg);
    allocInfo = ruInfo(cfg);
    ruSize = allocInfo.RUSizes(ruIdx);
elseif isa(cfg,'heTBSystemConfig')
    numHESIGB = 0;
    ruIdx = varargin{1};
    pktFormat = packetFormat(cfg);
    allocInfo = ruInfo(cfg);
    ruSize = allocInfo.RUSizes(ruIdx);
elseif isa(cfg,'wlanHERecoveryConfig')
    ruSize = cfg.RUSize;
    pktFormat = cfg.PacketFormat;
    if strcmp(pktFormat,'HE-MU')
        s = getSIGBLength(cfg);
        numHESIGB = s.NumSIGBSymbols;
    else
        numHESIGB = 0;
    end
else
    % SU, EXT SU, TB
    numHESIGB = 0;
    pktFormat = packetFormat(cfg);
    allocInfo = ruInfo(cfg);
    ruSize = allocInfo.RUSizes(ruIdx);
end

if strcmp(pktFormat,'HE-EXT-SU')
    numHESIGA = 4;
else % SU or MU
    numHESIGA = 2;
end

z = 2+numHESIGA+numHESIGB; % Pilot symbol offset
% Get the reference pilots for one space-time stream, pilot sequence same
% for all space-time streams
refPilots = wlan.internal.hePilots(ruSize,1,n,z); % Nsp-by-Nsym-by-1

% Average single-stream pilot estimates over symbols (2nd dimension)
avChanEstSSPilots = mean(chanEstSSPilots,2); % Nsp-by-1-by-Nrx

% Estimate channel at pilot location using least square estimates
chanEstPilotsLoc = bsxfun(@rdivide,x,refPilots); % Nsp-by-Nsym-by-Nrx

% Subtract the noisy least squares estimates of the channel at pilot symbol
% locations from the noise averaged single stream pilot symbol estimates of
% the channel
error = bsxfun(@minus,chanEstPilotsLoc,avChanEstSSPilots); % Nsp-by-Nsym-by-Nrx

% Get power of error and average over pilot symbols, subcarriers and
% receive antennas
useIdx = ~isnan(error); % NaNs may exist in 1xHELTF
nest = real(mean(error(useIdx).*conj(error(useIdx)),'all')); % For codegen

if nargout>1
    % Get power of channel estimate at pilot locations
   sigest =  real(mean(chanEstPilotsLoc(:).*conj(chanEstPilotsLoc(:)))); % For codegen
end

end