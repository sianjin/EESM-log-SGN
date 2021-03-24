function plotPERvsSNR(simParams,results)
% plotPERvsSNR Example helper function

%   Copyright 2019 The MathWorks, Inc.

% Get the unique channels, MIMO schemes and MCS simulated
uniqueChannels = unique([simParams.DelayProfile]);
uniqueNumTxRx = unique([[simParams.NumTransmitAntennas]; [simParams.NumReceiveAntennas]]','rows');
uniqueMCS = unique([simParams.MCS]);

for ichan = 1:numel(uniqueChannels)
    % Get the index of simulations results for the given delay profile
    channelIdx = [simParams.DelayProfile]==uniqueChannels(ichan);

    for inumtxrx = 1:size(uniqueNumTxRx,1)
        hi = gobjects(numel(uniqueMCS),1);
        ha = gobjects(numel(uniqueMCS),1);

        % Create a new figure for each channel and MIMO configuration
        fper = figure;
        axper = axes(fper); %#ok<LAXES>
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
            
            % If there are 
            snr = [simParams(resultIdx).SNR];
            resultsUse = [results{resultIdx}];

            per = [resultsUse.packetErrorRate];
            perAbs = [resultsUse.packetErrorRateAbs];
            
            cRow = mod(imcs-1,size(colors,1))+1; % Select color row
            hi(imcs) = semilogy(axper,snr,per,':*','Color',colors(cRow,:));
            hold on;
            ha(imcs) = semilogy(axper,snr,perAbs,'-o','Color',colors(cRow,:));
        end
        title(axper,['Abstracted vs Simulated PHY, ' num2str(uniqueNumTxRx(inumtxrx,1)) 'x' num2str(uniqueNumTxRx(inumtxrx,2)) ', ' char(uniqueChannels(ichan))]);
        grid(axper,'on')
        xlabel(axper,'SNR (dB)');
        ylabel(axper,'PER');
        
        % Discard any MCS plots which we don't have simulation results for
        huse = arrayfun(@(x) ~isa(x,'matlab.graphics.GraphicsPlaceholder'),ha);
        ha = ha(huse);
        hi = hi(huse);
        legend([ha; hi],[arrayfun(@(x) ['Abstracted, MCS ' num2str(x)],uniqueMCS(huse),'UniformOutput',false) arrayfun(@(x) ['Simulated, MCS ' num2str(x)],uniqueMCS(huse),'UniformOutput',false)], ...
            'location', 'best')
    end
end

end