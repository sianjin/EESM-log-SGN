%% MODEL-D, 1*1, MCS4, RU size 52, APEP length 1000
clear
load('rbirEffSnr_RU106_Model-D_8-by-2_MCS4.mat');

channelCoding = 'LDPC'; 
dataLength = 1000;  
format = 'HE_MU'; % hard code for MUConfig

abstraction = tgaxRBIRLinkPerformanceModel;
per = estimatePER(abstraction,snrEff,format,mcs,channelCoding,dataLength);
perAvg = nanmean(per,1); % Average PER, computed after removing all NaN values

heff = gobjects(1,1); 

% Create a new figure for each channel and MIMO configuration
fper = figure;
axeff = axes(fper); %#ok<LAXES>
colors = colormap('lines');
heff = semilogy(axeff,snrs,perAvg,'-o','Color',colors);

title(axeff,['Simulated link PER vs SNR, ' num2str(numTxRx(1)) 'x' num2str(numTxRx(2)) ', ' char(chan)]);
grid(axeff,'on')
xlabel(axeff,'SNR (dB)');
ylabel(axeff,'PER');

% Discard any MCS plots which we don't have simulation results for
huse = arrayfun(@(x) ~isa(x,'matlab.graphics.GraphicsPlaceholder'),heff);
heff = heff(huse);