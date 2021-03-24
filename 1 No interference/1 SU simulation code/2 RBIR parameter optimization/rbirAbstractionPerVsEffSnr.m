%% MODEL-D, 1*1, MCS4, RU size 52, APEP length 1000
clear
load('snrPer_CBW20_Model-D_1-by-1_MCS4.mat');

channelCoding = cfgHE.ChannelCoding; 
dataLength = cfgHE.APEPLength;  
format = 'HE_MU'; % hard code for MUConfig
BW = '20MHz';

abstraction = tgaxRBIRLinkPerformanceModel;

heff = gobjects(1,1); 

% Create a new figure for each channel and MIMO configuration
fper = figure;
axeff = axes(fper); %#ok<LAXES>
colors = colormap('lines');

resultIdx = logical(ones(1,size(results,2)));
resultsUse = [results{resultIdx}];
sinrStore = cat(3,resultsUse.sinrStore);
perStore = cat(1,resultsUse.perStore);
numSnr = sum(resultIdx);

% Initialize RBIR parameters
rbirParameters = [1 1]; % 1st parameter: alpha, 2nd parameter: beta
mse = @(rbirParameters)awgnPerSnrFittingMse(abstraction,sinrStore,perStore,format,mcs,channelCoding,dataLength,rbirParameters);
rbirParametersOpt = fminsearch(mse,rbirParameters);
[binsnr,binper,lut,snrEff] = awgnPerSnrFitting(abstraction,sinrStore,perStore,format,mcs,channelCoding,dataLength,rbirParametersOpt,numSnr);

% Store effective SNR vector under different SNR point
fname_I = sprintf('rbirEffSnr_%s_%s_%s-by-%s_MCS%s.mat',BW,char(chan),num2str(numTxRx(1)),num2str(numTxRx(2)),num2str(mcs));
save(fname_I,'snrEff','rbirParametersOpt','BW','mcs','numTxRx','chan','maxNumPackets','snrs')

cRow = 1; % Select color row
heff = semilogy(axeff,binsnr,binper,'x','Color',colors(cRow,:));
hold on;

href = semilogy(axeff,lut(:,1),lut(:,2),'-k');

title(axeff,['Simulated link PER vs effective SNR, ' num2str(numTxRx(1)) 'x' num2str(numTxRx(2)) ', ' char(chan)]);
grid(axeff,'on')
xlabel(axeff,'Effective SNR (dB)');
ylabel(axeff,'PER');

% Discard any MCS plots which we don't have simulation results for
huse = arrayfun(@(x) ~isa(x,'matlab.graphics.GraphicsPlaceholder'),heff);
heff = heff(huse);

legend([href; heff], [{'AWGN'}, arrayfun(@(x) ['MCS ' num2str(x)],mcs(huse),'UniformOutput',false)], ...
    'location', 'best')