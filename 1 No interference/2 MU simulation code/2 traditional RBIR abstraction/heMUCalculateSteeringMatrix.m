function steeringMatBF = heMUCalculateSteeringMatrix(steeringMatFB,cfg,cfgNDP,ruIdx)
%heMUCalculateSteeringMatrix Calculate beamforming steering matrix
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   STEERINGMATBF = heMUCalculateSteeringMatrix(STEERINGMATFB,CFG,CFGNDP,RUIDX) 
%   returns the steering matrix recommended to beamform an RU in a transmit
%   beamforming, or MU-MIMO configuration. ZF precoding is used.
%
%   STEERINGMATFB is a cell array containing the steering matrices fed-back
%   by each user in the RU to beamform.
%
%   CFG is the configuration of the HE-MU transmission and is a format
%   configuration object of type
%   <a href="matlab:help('wlanHEMUConfig')">wlanHEMUConfig</a>.
%   
%   CFGNDP is the configuration of the HE-NDP used to gather feedback and
%   is a format configuration object of type
%   <a href="matlab:help('wlanHESUConfig')">wlanHESUConfig</a>.
%
%   RUIDX is the RU index.

%   Copyright 2018 The MathWorks, Inc.

allocInfo = ruInfo(cfg);

% Indices of active subcarriers within the RU
ruOFDMInfo = wlanHEOFDMInfo('HE-Data',cfg,ruIdx);
ruInd = ruOFDMInfo.ActiveFrequencyIndices;

% Indices of active subcarriers in the NDP
ndpOFDMInfo = wlanHEOFDMInfo('HE-Data',cfgNDP);
trainingInd = ndpOFDMInfo.ActiveFrequencyIndices;

% Get the indices which overlap - use to extract from NDP
[~,scUseInd] = intersect(trainingInd,ruInd);

% Extract the RU of interest from the full bandwidth grid
numUsers = allocInfo.NumUsersPerRU(ruIdx);
steeringMatUse = cell(numUsers,1);

for i = 1:numUsers
    % Only take the RU subcarriers and space-time streams of
    % interest for the current RU and user
    userIdx = cfg.RU{ruIdx}.UserNumbers(i);
    numSTS = cfg.User{userIdx}.NumSpaceTimeStreams;
    numRx = size(steeringMatFB{userIdx},2);
    if numSTS>numRx
        error('The number of space-time streams (%d) exceeds the number of receive antennas (%d) for user %d',numSTS,numRx,userIdx);
    end
    steeringMatUse{i} = steeringMatFB{userIdx}(scUseInd,1:numSTS,:);
end

% Extract steering matrix for each RU
if numUsers>1
    steeringMatBF = muSteeringMatrixFromFeedback(steeringMatUse);
else
    steeringMatBF = steeringMatUse{1};
end

end

function steeringMatrix = muSteeringMatrixFromFeedback(mappingMatrix,varargin)
%   Q = muSteeringMatrixFromFeedback(QU) calculates the spatial mapping
%   matrix for a MU transmission using the ZF algorithm.
%
%   Q is an Nst-by-Nsts-by-Nt mapping matrix. Nst is the number of
%   subcarriers, Nsts is the total number of space-time streams, and Nt is
%   the number of transmit antennas.
%
%   QU is a cell array containing the individual mapping matrix for each
%   user. Each element of QU is sized Nst-by-Nstsu-by-Nt, where Nstsu is
%   the number of space-time streams for the individual user.
%
%   Q = muSteeringMatrixFromFeedback(QU,SNR) calculates the spatial mapping
%   matrix for a MU transmission using the MMSE algorithm given the SNR.
    if nargin>1
        precodingType = 'MMSE';
        snr = varargin{1};
    else
        precodingType = 'ZF';
    end

    numUsers = numel(mappingMatrix);

    % Get the number of STS per user
    numSTS = zeros(numUsers,1);
    for uIdx = 1:numUsers
        numSTS(uIdx) = size(mappingMatrix{uIdx},2);
    end
    numSTSTotal = sum(numSTS);

    % Pack the per user CSI into a matrix
    [numST,~,numTx] = size(mappingMatrix{1});        % Number of subcarriers
    steeringMatrix = zeros(numST,numTx,numSTSTotal); % Nst-by-Nt-by-Nsts

    for uIdx = 1:numUsers
        stsIdx = sum(numSTS(1:uIdx-1))+(1:numSTS(uIdx));
        steeringMatrix(:,:,stsIdx) = permute(mappingMatrix{uIdx},[1 3 2]); % Nst-by-Nt-by-Nsts
    end

    % Zero-forcing or MMSE precoding solution
    if strcmp(precodingType, 'ZF')
        delta = 0; % Zero-forcing
    else
        delta = (numTx/(10^(snr/10))) * eye(numTx); % MMSE
    end
    for i = 1:numST
        % Channel inversion precoding
        h = squeeze(steeringMatrix(i,:,:));
        steeringMatrix(i,:,:) = h/(h'*h + delta);
    end

    steeringMatrix = permute(steeringMatrix,[1 3 2]);
end