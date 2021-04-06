classdef tgaxEESMLinkPerformanceModel < handle
%tgaxLinkPerformanceModel Create a link performance model object
%   abstraction = tgaxLinkPerformanceModel returns a TGax link abstraction
%   model. This model is used to estimate the packet error rate for an
%   802.11ax single-user link assuming perfect synchronization.
%
%   tgaxEESMLinkPerformanceModel methods:
%
%   estimateLinkPerformance - returns the expected packet error rate given 
%                             the SINR per subcarrier.
%   effectiveSinrVec           - calculates the effective SINR given the SINR 
%                             per subcarrier.
%   estimatePER             - returns the estimate packet error rate given 
%                             the effective SINR.
%   selectAWGNLUT           - returns the appropriate AWGN lookup table.
%
%   % Example: estimate packet error rate for a link.
%
%   abstraction = tgaxLinkPerformanceModel;
%   sinr = 30*rand(224,2); % Number of subcarriers - by number of spatial streams
%   dataLen = 1e3; % Bytes
%   mcs = 3;
%   coding = 'BCC';
%   per = estimateLinkPerformance(abstraction,sinr,dataLen,mcs,coding)
%
%   See also calculateSINR.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

properties (Access=private)
    rbir; % RBIR LUT 
    awgn; % AWGN LUT
end

methods
    function obj = tgaxEESMLinkPerformanceModel(~)
        % Constructor, load lookup tables
        
        % RBIR table from 802.11-14/1450r0 - Box 0 Calibration Results
        lrbirlut = load('tgax_rbit.mat');
        obj.rbir = lrbirlut;
        
        % AWGN table generated with WLAN Toolbox
        awgnlut = load('tgax_awgn.mat');
        obj.awgn = awgnlut;
    end

    function [per,snreff] = estimateLinkPerformance(obj,sinr,varargin)
        % [PER,SNREFF] = estimateLinkPerformance(OBJ,SINR,DATALENGTH,FORMAT,MCS,CODING)
        % returns the estimated packet error rate and effective SINR.
        %
        % SINR is an array containing the SINR for each subcarrier,
        % symbol and spatial stream.
        %
        % DATALENGTH is the payload length in bytes.
        % 
        % FORMAT is one of 'NonHT','HTMixed','VHT','HE_SU','HE_EXT_SU'.
        %
        % MCS is the modulation and coding scheme index and must be
        % 0-9.
        %
        % CODING is the channel coding used and must be 'BCC' or 'LDPC'.
        %
        % [PER,SNREFF] = estimateLinkPerformance(OBJ,SINR,CFGSU) returns the
        % estimated packet error rate given the single-user format
        % configuration object CFGSU.
        %
        % [PER,SNREFF] = estimateLinkPerformance(OBJ,SINR,CFGMU,USERIDX)
        % returns the estimated packet error rate given the OFDMA
        % multi-user format configuration object CFGMU and user index
        % USERIDX.
        
        nargin = 4; % hard code for using MUConfig L2S mapping
        if nargin<6
            cfg = varargin{1};
            if nargin==4
                % [PER,SNREFF] = estimateLinkPerformance(OBJ,SINR,CFGMU,USERIDX)
                assert(isa(cfg,'wlanHEMUConfig'))
                % userIdx = varargin{2};
                userIdx = 1; % hard code for using MUConfig L2S mapping
                mcs = cfg.User{userIdx}.MCS;
                dataLength = cfg.User{userIdx}.APEPLength;
                coding = cfg.User{userIdx}.ChannelCoding;
                format = 'HE_MU';
            else
                % [PER,SNREFF] = estimateLinkPerformance(OBJ,SINR,CFGSU)
                mcs = cfg.MCS;

                % Get the channel coding and data long from the
                % configuration object
                switch class(cfg)
                    case {'wlanHESUConfig','wlanVHTConfig'}
                        dataLength = cfg.APEPLength;
                        coding = cfg.ChannelCoding;
                        format = 'HE_SU'; % Use HE-SU even if VHT as same MCS indices
                    case 'wlanHTConfig'
                        dataLength = cfg.PSDULength;
                        coding = cfg.ChannelCoding;
                        format = 'HTMixed';
                    case 'wlanNonHTConfig'
                        dataLength = cfg.PSDULength;
                        coding = 'BCC';
                        format = 'NonHT';
                    otherwise
                        error('Unexpected object');
                end
            end
        else
            % [PER,SNREFF] = estimateLinkPerformance(OBJ,SINR,DATALENGTH,FORMAT,MCS,CODING)
            narginchk(6,6)
            dataLength = varargin{1};
            format = varargin{2};
            mcs = varargin{3};
            coding = varargin{4};
        end

        % Calculate the effective SNR
        snreff = effectiveSINR(obj,sinr,format,mcs);

        % Estimate the packet error rate for the given SNR and
        % configuration
        per = estimatePER(obj,snreff,format,mcs,coding,dataLength);
    end

    function snrEff = effectiveSinrVec(obj,sinrStore,beta)
        % SNREFF = effectiveSINR(SINRSTORE,BETA) returns
        % the effective SNR vector.
        %
        % SINRSTORE is the SINR per subcarrier, spatial stream, and symbol.
        
        numPkts = size(sinrStore,3);
        snrEff = nan(numPkts,1);
        
        for pktIdx = 1: numPkts 
            sinr = sinrStore(:,:,pktIdx);
            sinr = 10.^(sinr/10);
            expAvg = mean(exp(-sinr/beta),'all');
            snreff = -beta*log(expAvg);
            snrEff(pktIdx) = 10*log10(snreff);
        end
        
    end

    function [per,perPL0,L0,lut] = estimatePER(obj,snreff,format,mcs,coding,dataLength)
        % [PER,PERPL0,L0,LUT] = estimatePER(OBJ,SNREFF,FORMAT,MCS,CODING,DATALENGTH)
        % returns the packet error rate PER, for the reference data length,
        % PERPL0, the reference data length L0, and the selected AWGN
        % lookup table, LUT.
        %
        % SNREFF is the effective SNR.
        %
        % FORMAT is one of 'NonHT','HTMixed','VHT','HE_SU','HE_EXT_SU'.
        %
        % MCS is the modulation and coding scheme index for the specified
        % format.
        %
        % CODING is either 'BCC' or 'LDPC'.
        %
        % DATALENGTH is the PSDU length in bytes.

        % Select AWGN lookup table
        [lut,L0] = selectAWGNLUT(obj,format,mcs,coding,dataLength);

        % Given the LUT and effective SNR, interpolate and extrapolate
        % to get the SNR
        perPL0 = interpolatePER(snreff,lut(:,:,1));

        % Calculate final PER, adjusting for packet size
        per = 1-(1-perPL0).^(dataLength/L0(1));
    end
    
    function [lut,L0] = selectAWGNLUT(obj,format,mcs,coding,dataLength)
        % [LUT,L0] = selectAWGNLUT(OBJ,FORMAT,MCS,CODING,DATALENGTH)
        % returns the appropriate AWGN lookup table, LUT, and reference
        % data length L0.
        %
        % LUT is a matrix containing the lookup table. Each row is of the
        % form [SNR PER].
        %
        % FORMAT is one of 'NonHT','HTMixed','VHT','HE_SU','HE_EXT_SU'.
        %
        % MCS is the HE modulation and coding scheme index and must be
        % 0-9.
        %
        % CODING is either 'BCC' or 'LDPC'.
        %
        % DATALENGTH is the PSDU length in bytes.

        [modulation, rate] = mcs2rate(format,mcs);
        switch coding
            case 'BCC'
                % As per TGax evaluation methodology, if the number of
                % bytes is <400 then use 32 byte BCC LUT, otherwise use
                % 1458 byte LUT.
                idx = all(bsxfun(@eq,obj.awgn.perTable_BCC_MCS,[modulation rate]),2);
                assert(any(idx),'Format and MCS is not valid with BCC coding')
                if dataLength < 400
                    lut = obj.awgn.perTable_BCC_32(:,:,idx);
                    L0 = obj.awgn.L0_BCC_32(idx);
                else
                    lut = obj.awgn.perTable_BCC_1458(:,:,idx);
                    L0 = obj.awgn.L0_BCC_1458(idx);
                end
            otherwise
                assert(strcmp(coding,'LDPC'))
                idx = all(bsxfun(@eq,obj.awgn.perTable_LDPC_MCS,[modulation rate]),2);
                assert(any(idx),'Format and MCS is not valid with LDPC coding')
                lut = obj.awgn.perTable_LDPC(:,:,idx);
                L0 = obj.awgn.L0_LDPC(idx);
        end
    end
end
end

function per = interpolatePER(snr,lut)
    % Lookup packet error rate, linear interpolation over logarithm LUT is
    % [SNR PER] matrix where each row is an entry
    
    % If there are multiple SNRs with PER=0 then ignore entries after the
    % first
    numSNRs = size(lut,1);
    idx = (1:numSNRs)';
    iszero = idx(lut(:,2)==0);
    if isempty(iszero)
        firstzero = numSNRs;
    else
        firstzero = iszero(1);
    end
    snrLUT = lut(1:firstzero,1);
    perLUT = lut(1:firstzero,2);
    
    per = 10.^(interp1(snrLUT,log10(perLUT),snr,'linear','extrap'));
    % If PER is less than 1 in a million then there are no errors.
    % This avoids the edge conditions when the effective SNR saturates due
    % to RBIR lookup causing a PER error floor (albeit very low) for high
    % SNRs.
    per(per<1e-6) = 0;
    per(per>1) = 1;
end

function [modscheme,rate] = mcs2rate(format,mcs)
    % Return the modulation scheme and rate given a format and MCS index
    switch format
        case 'NonHT'
            [modscheme,rate] = nonhtmcs2rate(mcs);
        case {'VHT','HE_SU','HE_EXT_SU','HE_MU'}
            [modscheme,rate] = hemcs2rate(mcs);
        otherwise
            assert(strcmp(format,'HTMixed'),'Invalid format')
            [modscheme,rate] = htmcs2rate(mcs);
    end
end

function [modscheme,rate] = hemcs2rate(mcs)
    % HE (and VHT) modulation and coding scheme from MCS
    switch mcs
        case 0
            modscheme = 2;
            rate = 1/2;
        case 1
            modscheme = 4;
            rate = 1/2;
        case 2
            modscheme = 4;
            rate = 3/4;
        case 3
            modscheme = 16;
            rate = 1/2;
        case 4
            modscheme = 16;
            rate = 3/4;
        case 5
            modscheme = 64;
            rate = 2/3;
        case 6
            modscheme = 64;
            rate = 3/4;
        case 7
            modscheme = 64;
            rate = 5/6;
        case 8
            modscheme = 256;
            rate = 3/4;
        case 9
            modscheme = 256;
            rate = 5/6;
        case 10
            modscheme = 1024;
            rate = 3/4;
        otherwise % 11
            assert(mcs==11);
            modscheme = 1024;
            rate = 5/6;
    end
end

function [modscheme,rate] = htmcs2rate(mcs)
    % HT modulation and coding scheme from MCS
    % HT MCS wraps around every 8 as it climbs spatial-streams
    [modscheme,rate] = hemcs2rate(mod(mcs,8));
end

function [modscheme,rate] = nonhtmcs2rate(mcs)
    % Non-HT modulation and coding scheme from MCS
    switch mcs
        case 0
            modscheme = 2;
            rate = 1/2;
        case 1
            modscheme = 2;
            rate = 3/4;
        case 2
            modscheme = 4;
            rate = 1/2;
        case 3
            modscheme = 4;
            rate = 3/4;
        case 4
            modscheme = 16;
            rate = 1/2;
        case 5
            modscheme = 16;
            rate = 3/4;
        case 6
            modscheme = 64;
            rate = 2/3;
        otherwise % 7
            assert(mcs==7,'Invalid MCS for Non-HT');
            modscheme = 64;
            rate = 3/4;           
    end
end