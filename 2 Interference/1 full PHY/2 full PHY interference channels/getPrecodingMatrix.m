function W = getPrecodingMatrix(cfg,varargin)
%getPrecodingMatrix(CFGHE) return the precoding matrix
%
%  W = getPrecodingMatrix(CFGHE) returns the precoding matrix W given the
%  format configuration object CFG.
%
%   W is a Nst-by-Nsts-by-Ntx precoding matrix, where Nst is the number of
%   active subcarriers, Nsts is the number of space-time streams, and Ntx
%   is the number of transmit antennas.
%
%   CFG is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, 
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, <a
%   href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>, <a
%   href="matlab:help('wlanS1GConfig')">wlanS1GConfig</a>, <a
%   href="matlab:help('wlanDMGConfig')">wlanDMGConfig</a>, or 
%   <a href="matlab:help('wlanHESUConfig')">wlanHESUConfig</a>.
%
%   W = getPrecodingMatrix(CFGHE,FIELD) returns the precoding used for the
%   field. FIELD can be 'data' or 'preamble'.
%
%   W = getPrecodingMatrix(CFGHEMU,RUIDX,FIELD) returns the precoding for
%   the RU specified by RUIDX.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

switch nargin
    case 1
        field = 'data';
    case 2
        if isnumeric(varargin{1})
            ruidx = varargin{1};
            field = 'data';
        else
            field = varargin{1};
            validatestring(field,{'data','preamble'},mfilename);
        end
    case 3
        ruidx = varargin{1};
        field = varargin{2};
        validatestring(field,{'data','preamble'},mfilename);
end

Wcs = getCyclicShiftMatrix(cfg,varargin{:});

if strcmp(field,'data') && ~isa(cfg,'wlanNonHTConfig')
    % Spatial mapping only relevant for:
    % * Data field, as not supporting BeamChange=false for now
    % * Configurations which perform spatial mapping 
    Wsm = getSpatialMappingMatrix(cfg,ruidx);
	W = bsxfun(@times,Wsm,Wcs); % Nst-by-Nsts-by-Ntx
else
    W = Wcs;
end

end

function csd = getCyclicShiftMatrix(cfg,varargin)
    % CSD = getCyclicShiftMatrix(CFG) returns a Nst-by-Nsts-by-1/Ntx matrix
    % containing the cyclic shift applied to each subcarrier and space-time
    % stream. Nst is the number of active subcarriers and Nsts is the
    % number of space-time streams. If the cyclic shift applied to each
    % transmitter is the same the size of third dimension returned is 1.
    %
    % CSD = getCyclicShiftMatrix(CFG,FIELD) the cyclic shift applied to the
    % FIELD. FIELD can be 'preamble' or 'data'.
    %
    % CSD = getCyclicShiftMatrix(CFGMU,RUIDX,FIELD) the cyclic shift
    % applied to an RU with index RUIDX.
    
    %#codegen
    
    ruIndx = 1;
    field = 'data';
    if nargin>1
        if isnumeric(varargin{1})
            ruIndx = varargin{1};
            if nargin>2
                field = varargin{2};
            end
        else
            field = varargin{1};
        end
    end

   info = ofdmInfo(cfg,field,ruIndx);
   cbw = wlan.internal.cbwStr2Num(cfg.ChannelBandwidth);

   % Get the cyclic shift either per space-time streams of per-transmit
   % antenna.
   % Note for the Non-HT format we always return the cyclic shift per
   % transmit antenna for the data and preamble
   switch field
       case 'data'
            switch class(cfg)
                case {'wlanHEMUConfig'}
                    allocInfo = ruInfo(cfg);
                    numSTS = allocInfo.NumSpaceTimeStreamsPerRU(ruIndx);
                    csh = wlan.internal.getCyclicShiftVal('VHT',allocInfo.NumSpaceTimeStreamsPerRU(ruIndx),cbw); % Same CSD for HE,VHT and HT
                    % Create 'mock' channel estimate to apply cyclic shift to.
                    in = ones(info.NumTones,numSTS);
                case {'wlanHESUConfig','wlanVHTConfig','wlanHTConfig'}
                    csh = wlan.internal.getCyclicShiftVal('VHT',cfg.NumSpaceTimeStreams,cbw); % Same CSD for HE,VHT and HT
                    % Create 'mock' channel estimate to apply cyclic shift to.
                    in = ones(info.NumTones,cfg.NumSpaceTimeStreams);
                case 'wlanNonHTConfig'
                    csh = wlan.internal.getCyclicShiftVal('OFDM',cfg.NumTransmitAntennas,cbw);
                    % Create 'mock' channel estimate to apply cyclic shift to.
                    in = ones(info.NumTones,cfg.NumTransmitAntennas);
            end
       otherwise % 'preamble'
           csh = wlan.internal.getCyclicShiftVal('OFDM',cfg.NumTransmitAntennas,cbw);
            % Create 'mock' channel estimate to apply cyclic shift to.
            in = ones(info.NumTones,cfg.NumTransmitAntennas);
   end
   
   % Apply per-sts or per-tx cyclic shift
   k = info.ActiveFrequencyIndices;
   csdTmp = wlan.internal.cyclicShiftChannelEstimate(in,csh,info.FFTLength,k);
   
   if strcmp(field,'preamble') || isa(cfg,'wlanNonHTConfig')
       csd = permute(csdTmp,[1 3 2]); % CSD applied over second dimension so permute to third dimension to represent transmit antennas
   else
       csd = csdTmp;
   end

end

function Q = getSpatialMappingMatrix(cfg,varargin)
%getSpatialMappingMatrix Returns spatial mapping matrix used.
%   Q = getSpatialMappingMatrix(CFG) returns the spatial mapping matrix
%   used for each occupied subcarrier in the data portion.
%
%   Q is Nst-by-Nsts-by-Ntx where Nst is the number of occupied
%   subcarriers, Nsts is the number of space-time streams, and Ntx is the
%   number of transmit antennas.
%
%   CFG is a format configuration object of type wlanHESUConfig,
%   wlanVHTConfig, or wlanHTConfig.
%
%   Q = getSpatialMappingMatrix(CFGMU,RUIDX) returns the spatial mapping
%   matrix for an HE-MU configuration. CFGMU is a wlanHEMUConfig object and
%   RUIDX is the index of the RU of interest.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

if isa(cfg,'wlanHEMUConfig')
    allocInfo = ruInfo(cfg);
    ruIdx = varargin{1};
    assert(ruIdx>0)
    numSTS = allocInfo.NumSpaceTimeStreamsPerRU(ruIdx);
    mappingType = cfg.RU{ruIdx}.SpatialMapping;
    mappingMatrix = cfg.RU{ruIdx}.SpatialMappingMatrix;
else
    ruIdx = 1; % Not used
    numSTS = sum(cfg.NumSpaceTimeStreams); % For VHT might be a vector
    mappingType = cfg.SpatialMapping;
    mappingMatrix = cfg.SpatialMappingMatrix;
end
numTx = cfg.NumTransmitAntennas;

info = ofdmInfo(cfg,'data',ruIdx);
Nst = info.NumTones;

switch mappingType
  case 'Direct'
    Q = repmat(permute(eye(numSTS,numTx),[3 1 2]),Nst,1,1);
  case 'Hadamard'
    hQ = hadamard(8);
    normhQ = hQ(1:numSTS,1:numTx)/sqrt(numTx);
    Q = repmat(permute(normhQ,[3 1 2]),Nst,1,1);
  case 'Fourier'
    [g1, g2] = meshgrid(0:numTx-1, 0:numSTS-1);
    normQ = exp(-1i*2*pi.*g1.*g2/numTx)/sqrt(numTx);
    Q = repmat(permute(normQ,[3 1 2]),Nst,1,1);
  otherwise  % case 'Custom'
    if size(mappingMatrix, 1) <= 8
        [numSTS,numTx] = size(mappingMatrix); % MappingMatrix is Nsts-by-Ntx
        Q = repmat(permute(normalize(mappingMatrix(1:numSTS, 1:numTx),numSTS),[3 1 2]),Nst,1,1);
    else
        Q = mappingMatrix;
        [Nst,numSTS,numTx] = size(Q); % MappingMatrix is Nst-by-Nsts-by-Ntx
        Qp = permute(Q,[2 3 1]);
        Qn = coder.nullcopy(complex(zeros(numSTS,numTx,Nst)));
        for i=1:Nst
            Qn(:,:,i) = normalize(Qp(:,:,i),numSTS); % Normalize mapping matrix. SHOULD WE DO THIS INTERNALLY?
        end
        Q = permute(Qn,[3 1 2]);
    end    
end

end

function Q = normalize(Q,numSTS)
    Q = Q * sqrt(numSTS)/norm(Q,'fro'); % Normalize mapping matrix
end

function info = ofdmInfo(cfg,field,varargin)

   if nargin>2
       ruIndx = varargin{1};
   end

   switch field
       case 'data'
           % Get OFDM info for data fields of formats
            switch class(cfg)
                case 'wlanHEMUConfig'
                    info = wlanHEOFDMInfo('HE-Data',cfg,ruIndx);
                case 'wlanHESUConfig'
                    info = wlanHEOFDMInfo('HE-Data',cfg);
                case 'wlanVHTConfig'
                    info = wlanVHTOFDMInfo('VHT-Data',cfg);
                case 'wlanHTConfig'
                    info = wlanHTOFDMInfo('HT-Data',cfg);
                case 'wlanNonHTConfig'
                    info = wlanNonHTOFDMInfo('NonHT-Data');
                otherwise
                    error('Unexpected format');
            end
       case 'preamble'
           % Get OFDM info for preamble fields of formats
            switch class(cfg)
                case 'wlanHEMUConfig'
                    info = wlanHEOFDMInfo('HE-SIG-A',cfg);
                case 'wlanHESUConfig'
                    info = wlanHEOFDMInfo('HE-SIG-A',cfg);
                case 'wlanVHTConfig'
                    info = wlanVHTOFDMInfo('VHT-SIG-A',cfg);
                case 'wlanHTConfig'
                    info = wlanHTOFDMInfo('HT-SIG',cfg);
                case 'wlanNonHTConfig'
                    info = wlanNonHTOFDMInfo('NonHT-Data');
                otherwise
                    error('Unexpected format');
            end
       otherwise
           error('Unexpected field')
   end

end