function simParams = getBox0SimParams(chans,numTxRx,mcs,cfgHE,maxNumErrors,maxNumPackets,userIdx,ruIdx,Nsts,intPathlossdB)
% getBox0SimParams Example helper function

%   Copyright 2019 The MathWorks, Inc.

% These arrays define the value and order SNRs are defined
channelConfigs = ["Model-B","Model-D"];
anteannaSNRConfigs = [1 1; 4 1; 4 2; 8 2];

snr = {
    % Model-B
    [ ...
        {... % 1x1
    [-3:4:9,11], ...   % MCS 0
    [1:4:13], ...     % MCS 1
    [2:4:18], ...     % MCS 2
    [7:4:19,21], ...     % MCS 3
    [9:4:25], ...     % MCS 4
    [13:4:29], ...    % MCS 5
    [14:4:30], ...    % MCS 6
    [16:4:32,34], ...    % MCS 7
    [18:4:34,36] ... % MCS 8
    [18:4:38] ...     % MCS 9
    21.5:4:37.5 ... % MCS 10
    22:4:38 ...     % MCS 11
    }; ...
    	{... % 4x1
    1:3:17, ...  % MCS 0
    5:3:23, ...  % MCS 1
    9:4:30, ...  % MCS 2
    12:4:36, ... % MCS 3
    16:4:40, ... % MCS 4
    20:4:43, ... % MCS 5
    22:4:46, ... % MCS 6
    24:4:48, ... % MCS 7
    26:4:50 ...  % MCS 8
    29:4:54 ...  % MCS 9
    21.5:4:37.5 ... % MCS 10
    22:4:38 ...     % MCS 11
    }; ...
        {... % 4x2
    1:3:17, ...  % MCS 0
    5:3:23, ...  % MCS 1
    9:4:30, ...  % MCS 2
    12:4:36, ... % MCS 3
    16:4:40, ... % MCS 4
    20:4:43, ... % MCS 5
    22:4:46, ... % MCS 6
    24:4:48, ... % MCS 7
    26:4:50 ...  % MCS 8
    29:4:54 ...  % MCS 9
    21.5:4:37.5 ... % MCS 10
    22:4:38 ...     % MCS 11
    }; ...
         {... % 8x2
    1:3:17, ...  % MCS 0
    5:3:23, ...  % MCS 1
    9:4:30, ...  % MCS 2
    12:4:36, ... % MCS 3
    16:4:40, ... % MCS 4
    20:4:43, ... % MCS 5
    22:4:46, ... % MCS 6
    24:4:48, ... % MCS 7
    26:4:50 ...  % MCS 8
    29:4:54 ...  % MCS 9
    21.5:4:37.5 ... % MCS 10
    22:4:38 ...     % MCS 11
    }; ...
    ];

% Model-D
    [ ...
        {... % 1x1
    -3:3:9, ...     % MCS 0
    0:4:12, ...     % MCS 1
    2:4:18, ...     % MCS 2
    8:4:20, ...     % MCS 3
    11:4:23, ...    % MCS 4
    [14:4:26,28], ...    % MCS 5
    16:4:28, ...    % MCS 6
    18:4:30, ...    % MCS 7
    21.5:4:37.5 ... % MCS 8
    20:4:36 ...     % MCS 9
    21.5:4:37.5 ... % MCS 10
    22:4:38 ...     % MCS 11
    }; ...
    	{... % 4x1
    -4:2:2, ...   % MCS 0
    5:2:15, ...   % MCS 1
    9:2:21, ...   % MCS 2
    12:2:24, ...  % MCS 3
    12:2:16, ...  % MCS 4
    20:3:32, ...  % MCS 5
    22:3:34, ...  % MCS 6
    24:3:39, ...  % MCS 7
    27:3:40 ...   % MCS 8
    29:3:44 ...   % MCS 9
    21.5:4:37.5 ... % MCS 10
    22:4:38 ...     % MCS 11
    }; ...
        {... % 4x2
    [11,13,14,15], ...   % MCS 0
    5:2:15, ...   % MCS 1
    9:2:21, ...   % MCS 2
    12:2:24, ...  % MCS 3
    [13 16 18 20], ...  % MCS 4
    20:3:32, ...  % MCS 5
    22:3:34, ...  % MCS 6
    24:3:39, ...  % MCS 7
    27:3:40 ...   % MCS 8
    29:3:44 ...   % MCS 9
    21.5:4:37.5 ... % MCS 10
    22:4:38 ...     % MCS 11
    }; ...
        {... % 8x2
    1:2:10, ...   % MCS 0
    5:2:15, ...   % MCS 1
    9:2:21, ...   % MCS 2
    12:2:24, ...  % MCS 3
    [10,14,17,19], ...  % MCS 4
    20:2:30, ...  % MCS 5
    22:3:34, ...  % MCS 6
    24:3:39, ...  % MCS 7
    27:3:40 ...   % MCS 8
    29:3:44 ...   % MCS 9
    21.5:4:37.5 ... % MCS 10
    22:4:38 ...     % MCS 11
    }; ...
    ] ...
    };

% Create channel configuration
tgaxChannel = wlanTGaxChannel;
tgaxChannel.DelayProfile = 'Model-D';
tgaxChannel.NumTransmitAntennas = cfgHE.NumTransmitAntennas;
tgaxChannel.NumReceiveAntennas = 1;
tgaxChannel.TransmitReceiveDistance = 15; % Distance in meters for NLOS
tgaxChannel.ChannelBandwidth = cfgHE.ChannelBandwidth;
tgaxChannel.LargeScaleFadingEffect = 'None';
fs = wlanSampleRate(cfgHE);
tgaxChannel.SampleRate = fs;
tgaxChannel.PathGainsOutputPort = true;
tgaxChannel.NormalizeChannelOutputs = false;

% Generate a structure array containing the simulation parameters,
% simParams. Each element contains the parameters for a simulation.
simParamsRef = struct('MCS',0,'SNR',0,'RandomSubstream',0,'Config',cfgHE, ...
    'MaxNumPackets',maxNumPackets,'MaxNumErrors',maxNumErrors, ...
    'NumTransmitAntennas',0,'NumReceiveAntennas',0,'DelayProfile',"Model-B",...
    'Channel',tgaxChannel,'userIdx',userIdx,'ruIdx',ruIdx,'intPathlossdB',intPathlossdB);
simParams = repmat(simParamsRef,0,0);
% There must be a SNR cell for each channel
assert(all(numel(channelConfigs)==numel(snr)))
% There must be a SNR cell element for each MIMO configuration
assert(all(size(anteannaSNRConfigs,1)==cellfun(@(x)size(x,1),snr)))
for ichan = 1:numel(chans)
    channelIdx = chans(ichan)==channelConfigs;
    for itxrx = 1:size(numTxRx,1)
        numTxRxIdx = all(numTxRx(itxrx,:)==anteannaSNRConfigs,2);
        for imcs = 1:numel(mcs)
            snrIdx = mcs(imcs)+1;
            for isnr = 1:numel([snr{channelIdx}{numTxRxIdx,snrIdx}])
                % Set simulation specific parameters
                sp = simParamsRef;
                sp.MCS = mcs(imcs);
                sp.NumTransmitAntennas = numTxRx(itxrx,1);
                sp.NumReceiveAntennas = numTxRx(itxrx,2);
                sp.DelayProfile = chans(ichan);

                % Set random substream for reproducible results
                sp.RandomSubstream = isnr;
                
                % Setup PHY configuration for MUConfig
                for userIdx = 1:numel(sp.Config.User)
                    sp.Config.User{userIdx}.MCS = mcs(imcs);
                    sp.Config.NumTransmitAntennas = numTxRx(itxrx,1);
                    sp.Config.User{userIdx}.NumSpaceTimeStreams = Nsts;   
                end
                for ruIdxIter = 1:numel(sp.Config.RU)
                    % Initialize spatial mapping to be Fourier
                    sp.Config.RU{ruIdxIter}.SpatialMapping = 'Fourier'; 
                end             

                % Configure channel model
                sp.Channel = clone(tgaxChannel);
                sp.Channel.DelayProfile = chans(ichan);
                sp.Channel.NumTransmitAntennas = numTxRx(itxrx,1);
                sp.Channel.NumReceiveAntennas = numTxRx(itxrx,2);

                % Lookup SNR to simulate
                sp.SNR = snr{channelIdx}{numTxRxIdx,snrIdx}(isnr);

                % Append to other tests
                simParams = [simParams sp]; %#ok<AGROW>
            end
        end
    end
end

end