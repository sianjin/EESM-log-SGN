function [y,csi] = heEqualizeCombine(x,chanEst,nVar,cfg,varargin)
%heEqualizeCombine HE MIMO channel equalization and STBC combining
%
%   [Y,CSI] = heEqualizeCombine(X,CHANEST,NOISEVAR,CFGHE) performs
%   minimum-mean-square-error (MMSE) frequency domain equalization and
%   optionally STBC combining using the signal input X, the channel
%   estimate, CHANEST, and noise variance, NVAR.
%
%   Y is an estimate of the transmitted frequency domain signal and is of
%   size Nsd-by-Nsym-by-Nsts, where Nsd represents the number of carriers
%   (frequency domain), Nsym represents the number of symbols (time
%   domain), and Nsts represents the number of space-time streams (spatial
%   domain). It is complex when either X or CHANEST is complex, or is real
%   otherwise.
%
%   CSI is a real matrix of size Nsd-by-Nsts containing the soft channel
%   state information.
%
%   X is a real or complex array containing the frequency domain signal to
%   equalize. It is of size Nsd-by-Nsym-by-Nr, where Nr represents the
%   number of receive antennas.
%
%   CHANEST is a real or complex array containing the channel estimates for
%   each carrier and symbol. It is of size Nsd-by-Nsts-by-Nr.
%
%   NVAR is a nonnegative scalar representing the noise variance.
%
%   CFGHE is the format configuration object of type <a href="matlab:help('wlanHESUConfig')">wlanHESUConfig</a>,
%   <a href="matlab:help('wlanHETBConfig')">wlanHETBConfig</a> or <a href="matlab:help('wlanHERecoveryConfig')">wlanHERecoveryConfig</a>.
%
%   [Y,CSI] = heEqualizeCombine(X,CHANEST,NVAR,CFGMU,USERIDX) performs
%   equalization of a HE multi user transmission. The user index is
%   required to extract the space-time streams of interest for the user.
%
%   CFGMU is the format configuration object of type <a href="matlab:help('wlanHEMUConfig')">wlanHEMUConfig</a> or,
%   <a href="matlab:help('heTBSystemConfig')">heTBSystemConfig</a>, which specifies the parameters for multi user
%   HE format and trigger-based (HT TB) system format configuration object
%   respectively.
%
%   USERIDX is the index of the user to decode within the RU.

%   Copyright 2017-2019 The MathWorks, Inc.

%#codegen

validateattributes(cfg,{'wlanHESUConfig','wlanHEMUConfig','wlanHETBConfig','heTBSystemConfig','wlanHERecoveryConfig'},{'scalar'},mfilename,'format configuration object');

% Defaults
userIdx = 1;
if isa(cfg,'wlanHEMUConfig')
    narginchk(5,5) % Require user index for multi-user reception
    if nargin>4
        userIdx = varargin{1};
    end
    % Get the indices of the space-time streams for this user
    allSTSIdx = wlan.internal.heSpaceTimeStreamIndices(cfg);
    stsIdx = allSTSIdx(1,userIdx):allSTSIdx(2,userIdx);
elseif isa(cfg,'wlanHERecoveryConfig') && strcmp(cfg.PacketFormat,'HE-MU')
    % Get the indices of the space-time streams for this user
    stsIdx = cfg.SpaceTimeStreamStartingIndex:(cfg.SpaceTimeStreamStartingIndex+cfg.NumSpaceTimeStreams-1);
else
    stsIdx = 1; % For codegen
end

if cfg.STBC 
    % Only SU, get num of SS from size of channel estimate
    nss = size(chanEst,2)/2;
    [y,csi] = wlan.internal.wlanSTBCCombine(x,chanEst,nss,'MMSE',nVar);
else
    % Equalize 
    [y,csi] = helperSymbolEqualize(x,chanEst,nVar);
    if isa(cfg,'wlanHEMUConfig') || (isa(cfg,'wlanHERecoveryConfig') && strcmp(cfg.PacketFormat,'HE-MU'))
        % Extract used STS
        y = y(:,:,stsIdx);
        csi = csi(:,stsIdx);
    end
end
end