function plotPERvsSNR(simParams,results,numTxRx,chan)

snrCell = {simParams.SNR};
snr =[];
per = [];
for snrIdx = 1:numel(simParams)
    per = [per, results{snrIdx}.packetErrorRate];
    snr= [snr, snrCell{snrIdx}];
end

heff = gobjects(1,1);

fper = figure;
axeff = axes(fper); %#ok<LAXES>
colors = colormap('lines');

cRow = 1; % Select color row
href = semilogy(axeff,snr,per,'-k');

title(axeff,['Simulated link PER vs RX SNR, ' num2str(numTxRx(1)) 'x' num2str(numTxRx(2)) ', ' char(chan)]);
grid(axeff,'on')
xlabel(axeff,'Effective SNR (dB)');
ylabel(axeff,'PER');

% Discard any MCS plots which we don't have simulation results for
huse = arrayfun(@(x) ~isa(x,'matlab.graphics.GraphicsPlaceholder'),heff);
heff = heff(huse);

end