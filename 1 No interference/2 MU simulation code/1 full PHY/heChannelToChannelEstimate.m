function chanEst = heChannelToChannelEstimate(chan,cfg,varargin)
%heChannelToChannelEstimate Return the perfect channel estimate for HE fields
%   CHANEST = heChannelToChannelEstimate(CHAN,CFGSU) incorporates any
%   precoding applied, and the cyclic shifts per space-time stream to form
%   the perfect channel estimate, CHANEST, as estimated using the HE-LTF
%   fields.
%
%   CHANEST is a Nst-by-Nsym-by-Nsts-by-Nr array where Nst is the number of
%   active (occupied) subcarriers, Nsym is the number of symbols, Nsts is
%   the number of space-time streams, and Nr is the number of receive
%   antennas.
%
%   CHAN is the frequency domain channel response for active subcarriers
%   and is a Nst-by-Nsym-by-Nt-by-Nr array where Nt is the number of
%   transmit antennas.
%
%   CFGSU is a single-user format configuration object of type <a href="matlab:help('wlanHESUConfig')">wlanHESUConfig</a>.
%
%   CHANEST = heChannelToChannelEstimate(CHAN,CFGMU,RUOFINTEREST) returns
%   the perfect channel estimate for the RU of interest index RUOFINTEREST
%   for a multi-user configuration.
%
%   CFGMU is a single-user format configuration object of type <a href="matlab:help('wlanHEMUConfig')">wlanHEMUConfig</a>.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

W = getPrecodingMatrix(cfg,varargin{:}); % Nst-by-Nsts-by-Ntx

% Apply spatial mapping and cyclic shift
tmp = sum(bsxfun(@times,chan,permute(W,[1 4 3 5 2])),3); % Nst-by-Nsym-by-1-by-Nr-by-Nsts

% Permute to Nst-by-Nsym-by-Nsts-by-Nr
chanEst = permute(tmp,[1 2 5 4 3]);

% Scale the channel estimate based on the scaling applied at the
% transmitter and demodulator to match the practical channel estimator.
if nargin>2
    ruIdx = varargin{1};
    allocInfo = ruInfo(cfg);
    numSTS = allocInfo.NumSpaceTimeStreamsPerRU;
    alpha = allocInfo.PowerBoostFactorPerRU;
    ruSize = allocInfo.RUSizes;
    ruScalingFactor = sqrt(ruSize(ruIdx))*alpha(ruIdx)/sqrt(numSTS(ruIdx));
    allScalingFactor = 1/sqrt(sum(alpha.^2.*ruSize));
    chanEst = chanEst*allScalingFactor*ruScalingFactor;
else
    chanEst = chanEst/sqrt(sum(cfg.User{1}.NumSpaceTimeStreams));
end

end