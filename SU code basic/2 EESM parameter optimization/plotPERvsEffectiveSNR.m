function plotPERvsEffectiveSNR(simParams,results,varargin)
% plotPERvsEffectiveSNR Example helper function

%   Copyright 2019 The MathWorks, Inc.

if nargin>2
    abstraction = varargin{1};
else
    abstraction = tgaxLinkPerformanceModel;
end

% Get the unique channels, MIMO schemes and MCS simulated
uniqueChannels = unique([simParams.DelayProfile]);
uniqueNumTxRx = unique([[simParams.NumTransmitAntennas]; [simParams.NumReceiveAntennas]]','rows');
uniqueMCS = unique([simParams.MCS]);
switch class(simParams(1).Config)
    case 'wlanHESUConfig'
        format = 'HE_SU';
    case 'wlanVHTConfig'
        format = 'VHT';
    case 'wlanHTConfig'
        format = 'HTMixed';
    case 'wlanNonHTConfig'
        format = 'NonHT';
    case 'wlanHEMUConfig' % added for OFDMA L2S mapping
        format = 'HE_MU';
end

for ichan = 1:numel(uniqueChannels)
    % Get the index of simulations results for the given delay profile
    channelIdx = [simParams.DelayProfile]==uniqueChannels(ichan);

    for inumtxrx = 1:size(uniqueNumTxRx,1)
        heff = gobjects(numel(uniqueMCS),1);

        % Create a new figure for each channel and MIMO configuration
        fper = figure;
        axeff = axes(fper); %#ok<LAXES>
        colors = colormap('lines');

        % Get the index of simulation results for the given MIMO scheme
        ntxrxIdx = all([[simParams.NumTransmitAntennas]; [simParams.NumReceiveAntennas]]'==uniqueNumTxRx(inumtxrx,:),2)';
        
        for imcs = 1:numel(uniqueMCS)
            % Get the index for simulation results for the given MCS
            mcsIdx = [simParams.MCS]==uniqueMCS(imcs);

            % Extract results for this particular combination
            resultIdx = channelIdx&ntxrxIdx&mcsIdx;
            
            % If a result is empty then the simulation did not run, so
            % remove from result index
            for ir = 1:numel(resultIdx)
                resultIdx(ir) = resultIdx(ir)&~isempty(results{ir});
            end
            
            if ~any(resultIdx)
                % If there are no results to display then try next one
                continue
            end
            
            resultsUse = [results{resultIdx}];
            snrEff = [resultsUse.snreffStore];
            perStore = [resultsUse.perStore];
                       
            % Bin effective SNRs (0.25dB steps as per TGax EM) and calculate the PER in each bin
            edges = floor(min(snrEff(:))):0.25:ceil(max(snrEff(:)));
            BINS = discretize(snrEff(:),edges);
            binper = zeros(numel(edges)-1,1);
            for i = 1:numel(edges)-1
                binper(i) = sum(perStore(BINS==i))/sum(BINS==i);
            end

            cRow = mod(imcs-1,size(colors,1))+1; % Select color row
            heff(imcs) = semilogy(axeff,edges(1:end-1)+(edges(2)-edges(1))/2,binper,'x','Color',colors(cRow,:));
            hold on;
            
            % Plot AWGN reference curve for MCS, assume channel coding and
            % data length the same for all simulations
            % channelCoding = [simParams(1).Config.ChannelCoding];
            channelCoding = 'LDPC'; % hard code for MUConfig
            dataLength = [simParams(1).Config.User{1}.APEPLength];  % hard code for MUConfig
            lut = selectAWGNLUT(abstraction,format,uniqueMCS(imcs),channelCoding,dataLength);
            href = semilogy(axeff,lut(:,1),lut(:,2),'-k'); 
        end
        
        title(axeff,['Simulated link PER vs effective SNR, ' num2str(uniqueNumTxRx(inumtxrx,1)) 'x' num2str(uniqueNumTxRx(inumtxrx,2)) ', ' char(uniqueChannels(ichan))]);
        grid(axeff,'on')
        xlabel(axeff,'Effective SNR (dB)');
        ylabel(axeff,'PER');
        
        % Discard any MCS plots which we don't have simulation results for
        huse = arrayfun(@(x) ~isa(x,'matlab.graphics.GraphicsPlaceholder'),heff);
        heff = heff(huse);
        
        legend([href; heff], [{'AWGN'}, arrayfun(@(x) ['MCS ' num2str(x)],uniqueMCS(huse),'UniformOutput',false)], ...
            'location', 'best')
    end
end

end