function H = helperPerfectChannelEstimate(pathGains,pathFilters,varargin)
%helperPerfectChannelEstimate perfect channel estimation
%   H = helperPerfectChannelEstimate(PATHGAINS,PATHFILTERS,NFFT,NCP,ACTIVEFFTIND)
%   performs perfect channel estimation.
%
%   H is an array of size Nst-by-Nsym-by-Nt-by-Nr-by-Nl. Nst is the number
%   of active subcarriers. Nsym is the number of OFDM symbols. Nt is the
%   number of transmit antennas. Nr is the number of receive antennas. Nl
%   is the number of links.
%
%   PATHGAINS must be an array of size Ns-by-Np-by-Nt-by-Nr-by-Nl, where Ns
%   is the number of path gain samples and, Np is the number of paths. The
%   channel impulse response is averaged across all samples and summed
%   across all transmit antennas and receive antennas before timing
%   estimation.
%
%   PATHFILTERS must be a matrix of size Np-by-Nh where Nh is the number of
%   impulse response samples. The path filters is assumed to be the same
%   for all links.
%
%   NFFT is the FFT length.
%
%   NCP is the cyclic prefix length.
%
%   ACTIVEFFTIND is a vector of active subcarrier indices within the FFT
%   (range 1:NFFT).
%
%   H = helperPerfectChannelEstimate(...,OFFSET) performs perfect channel
%   estimation given a timing offset, OFFSET. If not provided the ideal
%   offset is calculated internally.
%
%   OFFSET is a vector of length Nl indicating estimated timing offset, an
%   integer number of samples relative to the first sample of the channel
%   impulse response reconstructed from PATHGAINS and PATHFILTERS.
%
%   H = helperPerfectChannelEstimate(...,'OFDMSymbolOffset',SYMOFFSET)
%   specifies the optional OFDM symbol sampling offset as a fraction of the
%   cyclic prefix length between 0 and 1, inclusive. When unspecified, a
%   value of 0 is used.

%   See also channelDelay.

%   Copyright 2019 The MathWorks, Inc. 

%#codegen

cpFraction = 0; % Default OFDM symbol offset

if nargin>5
    if isnumeric(varargin{4})
        % helperPerfectChannelEstimate(pathGains,pathFilters,FFTLen,CPLen,ActiveFFTInd,[Offset],[Name],[Value])
        offset = varargin{4};
        cpFraction = demodNVPairParse(varargin{5:end});
    else
        % helperPerfectChannelEstimate(pathGains,pathFilters,FFTLen,CPLen,ActiveFFTInd,[Name],[Value])
        cpFraction = demodNVPairParse(varargin{4:end});
    end
else
    % helperPerfectChannelEstimate(pathGains,pathFilters,FFTLen,CPLen,ActiveFFTInd)
    offset = channelDelay(pathGains,pathFilters);
end
ofdmInfo = struct;
ofdmInfo.FFTLength = varargin{1};
ofdmInfo.CPLength = varargin{2};
ofdmInfo.ActiveFFTIndices = varargin{3};
ofdmInfo.NumTones = numel(varargin{3});

validateInputs(pathGains,pathFilters);

sampleIndex = (1:(size(pathGains,1))).';

% Get number of channel impulse response samples 'Nh'
Nh = size(pathFilters,2);

% Get number of paths 'Np', number of transmit antennas 'Nt' and number
% of receive antennas 'Nr' in the path gains array. The pathGains are
% of size Ns-by-Np-by-Nt-by-Nr-by-Nl, where 'Ns' is the number of channel
% snapshots and Nl is the number of links.
[Ns,Np,Nt,Nr,Nl] = size(pathGains);

cpLength = ofdmInfo.CPLength(1);
fftSize = ofdmInfo.FFTLength;
activeFFTIndices = ofdmInfo.ActiveFFTIndices;
Nst = ofdmInfo.NumTones;

% Establish the starting and ending sample indices of each OFDM symbol
% across the total number of subframes, taking into consideration the
% initial slot number, and update the cyclic prefix lengths to span all
% subframes.
% Establish how many OFDM symbols 'L' are spanned by 'T' time samples.

% Return OFDM symbols for all symbols worth of data passed
symLength = (fftSize+cpLength);
L = ceil(Ns/symLength); % number of OFDM symbols
inc = 0:symLength:((L-1)*symLength);
symbolStarts = inc+cpLength;
symbolEnds = symbolStarts+fftSize;

% Ensure that total number of samples 'Ns' is at least symbol
Ns = max(Ns,symbolEnds(end));

% Path gains Ns-by-Np-by-Nt-by-Nr-by-Nl
% Path filters Nh-by-Np-by-Nl
H = complex(zeros(Nst,L,Nr,Nt,Nl,'like',pathGains));

for nl = 1:Nl
    symbolStartIdx = symbolStarts + offset(nl);
    idx = sum(bsxfun(@(x,y)x>=y,symbolStartIdx,sampleIndex),1);

    % Prepare the path gains matrix by indexing using 'idx' to select a
    % first dimension element for each OFDM symbol start, and permute to
    % put the multipath components in the first dimension and switch the
    % antenna dimensions. The pathGains are now of size Np-by-Nr-by-Nt
    pathGainsLink = pathGains(idx,:,:,:,nl);
    pathGainsLink = permute(pathGainsLink,[2 4 3 1]);

    % Create channel impulse response array 'h' for each impulse response
    % sample, receive antenna, transmit antenna and OFDM symbol
    h = zeros(Nh,Nr,Nt,L,'like',pathGainsLink);

    % For each path, add its contribution to the channel impulse response
    % across all transmit antennas, receive antennas and OFDM symbols
    for np = 1:Np
        h = h + bsxfun(@times,pathFilters(np,:).',pathGainsLink(np,:,:,:));
    end

    % Create the empty received waveform (for each transmit antenna)
    rxWave = zeros([Ns Nr Nt],'like',pathGainsLink);

    % For each OFDM symbol, add the corresponding impulse response samples
    % across all transmit antennas and receive antennas to the received
    % waveform. Note that the impulse responses are positioned according to
    % the timing offset 'offset' and the channel filter delay so that
    % channel estimate produced is as similar as possible to that produced
    % for a filtered waveform (without incurring the time cost of the full
    % filtering)
    for l = 1:L
        tl = symbolStarts(l) - offset(nl) + (1:Nh);
        rxWave(tl,:,:) = rxWave(tl,:,:) + h(:,:,:,l);
    end

    % Remove any samples from the end of the received waveforms that
    % correspond to incomplete OFDM symbols
    rxWave = rxWave(1:symbolEnds(L),:,:);

    symOffset = fix(cpLength * cpFraction);

    tstart = min(symbolStarts-cpLength)+1;
    for nt = 1:Nt
        rxGrid = ofdmdemod(rxWave(tstart:end,:,nt),fftSize,cpLength,symOffset);
        H(:,:,:,nt,nl) = rxGrid(activeFFTIndices,1:L,1:Nr); % Nst-by-Nsym-by-Nr (indexing for codegen)        
    end
end

H = permute(H,[1 2 4 3 5]); % Nst-by-Nsym-by-Nt-by-Nr-by-Nl

end

function validateInputs(pathGains,pathFilters)
% Check inputs
    
    % Validate channel path gains
    assert(~(ndims(pathGains)>5),'The number of dimensions in the path gains %d must be less than or equal to 5',ndims(pathGains));
    
    % Validate path filters impulse response
    assert(~(size(pathGains,2)~=size(pathFilters,1)),'The number of paths (2nd dimension size) in the path gains %d must equal the number of paths (1st dimension size) in the path filters %d.',size(pathGains,2),size(pathFilters,2));
    
end

function symOffset = demodNVPairParse(varargin)
%demodNVPairParse Parse optional name-value pairs for demodulation
% 
%   SYMOFFSET = demodNVPairParse(IN) parses the cell array IN to determine
%   the optional OFDM symbol offset as a name-value pair.

    % Default
    symOffset = 0;

    if nargin==0
        % No need for NV parsing
        return
    end

    % Validate each P-V pair
    if isempty(coder.target) % Simulation path
        p = inputParser;
        % Set defaults for the optional arguments
        addParameter(p,'OFDMSymbolOffset',symOffset);
        parse(p,varargin{:}); % Parse inputs
        res = p.Results; 
        symOffset = res.OFDMSymbolOffset;
    else % Codegen path
        if ~(nargin==2 && coder.internal.isConst(varargin{1}) && coder.internal.isTextRow(varargin{1}))
            % Error if name is not constant or NV pair npt provided
            coder.internal.error('wlan:shared:InvalidDemodNV');
            return
        end
        pvPairs = struct('OFDMSymbolOffset',uint32(0));
        % Select parsing options
        popts = struct('PartialMatching',true);
        % Parse inputs
        pStruct = coder.internal.parseParameterInputs(pvPairs,popts,varargin{:});
        % Get values for the P-V pair or set defaults for the optional arguments
        symOffset = coder.internal.getParameterValue(pStruct.OFDMSymbolOffset,symOffset,varargin{:});
    end
    validateattributes(symOffset,{'numeric'},{'scalar','>=',0,'<=',1},mfilename,'OFDM symbol offset');

end
