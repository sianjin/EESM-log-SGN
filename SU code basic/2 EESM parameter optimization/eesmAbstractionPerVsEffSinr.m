%% Introduction 
% This file is the main script to optimize the EESM parameter (beta)
% shown in Fig.6 of the IEEE TCOM paper:
% "Efficient PHY Layer Abstraction for Fast Simulations in Complex 
% System Environments"
%% Load data
clear all
% Copy the generated file from the full PHY simualtion
% in this folder and load it here
load('snrPer_CBW20_Model-D_2-by-4_MCS0.mat');

% Tuning parameter in this function: EESM parameter - beta
% Suggestion: the higher MCS, the larger initialized beta
% Initial values for reference: 
%   MCS: 0  -> beta : 1
%   MCS: 1  -> beta : 2
%   MCS: 2  -> beta : 1.5
%   MCS: 3  -> beta : 5
%   MCS: 4  -> beta : 7
%   MCS: 5  -> beta : 26
%   MCS: 6  -> beta : 33
%   MCS: 7  -> beta : 43
%   MCS: 8  -> beta : 111
%   MCS: 9  -> beta : 170
%   MCS: 10  -> beta : 410
%   MCS: 11  -> beta : 650
beta = 1;

% Load basic setup
channelCoding = cfgHE.ChannelCoding; 
dataLength = cfgHE.APEPLength;  
format = 'HE_MU'; % hard code for MUConfig
bandwidth = cfgHE.ChannelBandwidth;
abstraction = tgaxEESMLinkPerformanceModel;
 
% Load post-MIMO processing SINR matrix and error state
resultIdx = logical(ones(1,size(results,2)));
resultsUse = [results{resultIdx}];
sinrStore = cat(3,resultsUse.sinrStore);
perStore = cat(1,resultsUse.perStore);
numSnr = sum(resultIdx);

%% Optimize EESM parameter beta
mse = @(beta)awgnPerSnrFittingMse(abstraction,sinrStore,perStore,format,mcs,channelCoding,dataLength,beta);
betaOpt = fminsearch(mse,beta); % Optimized EESM parameter
[binsnr,binper,lut,snrEff] = awgnPerSnrFitting(abstraction,sinrStore,perStore,format,mcs,channelCoding,dataLength,betaOpt,numSnr);

%% Store effective SNR vector under different SNR point
fname_I = sprintf('eesmEffSinr_%s_%s_%s-by-%s_MCS%s.mat',bandwidth,char(chan),num2str(numTxRx(1)),num2str(numTxRx(2)),num2str(mcs));
save(fname_I,'snrEff','betaOpt','bandwidth','mcs','numTxRx','chan','maxNumPackets','snrs')

%% Plot EESM L2S mapping results
heff = gobjects(1,1);
fper = figure;
axeff = axes(fper); %#ok<LAXES>
colors = colormap('lines');
cRow = 1; % Select color row
heff = semilogy(axeff,binsnr,binper,'x','Color',colors(cRow,:));
hold on;
href = semilogy(axeff,lut(:,1),lut(:,2),'-k');
title(axeff,['Simulated Instantaneous PER vs Effective SINR, ' num2str(numTxRx(1)) 'x' num2str(numTxRx(2)) ', ' char(chan)]);
grid(axeff,'on')
xlabel(axeff,'Effective SINR (dB)');
ylabel(axeff,'Instataneous PER');
% Discard any MCS plots which we don't have simulation results 
huse = arrayfun(@(x) ~isa(x,'matlab.graphics.GraphicsPlaceholder'),heff);
heff = heff(huse);
legend([href; heff], [{'AWGN-SISO'}, {'EESM'},'location', 'best'])