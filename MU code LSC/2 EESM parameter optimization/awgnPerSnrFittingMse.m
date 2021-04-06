function mse = awgnPerSnrFittingMse(abstraction,sinrStore,perStore,format,mcs,channelCoding,dataLength,beta)

snrEff = effectiveSinrVec(abstraction,sinrStore,beta);

finiteSnrEffIdx = isfinite(snrEff);
snrEff = snrEff(finiteSnrEffIdx);
perStore = perStore(finiteSnrEffIdx);

% Bin effective SNRs (0.25dB steps as per TGax EM) and calculate the PER in each bin
edges = floor(min(snrEff(:))):0.25:ceil(max(snrEff(:)));
binsnr = edges(1:end-1)+(edges(2)-edges(1))/2;
BINS = discretize(snrEff(:),edges);
binper = zeros(numel(edges)-1,1);
for i = 1:numel(edges)-1
    binper(i) = sum(perStore(BINS==i))/sum(BINS==i);
end

% Plot AWGN reference curve for MCS, assume channel coding and
% data length the same for all simulations
lut = selectAWGNLUT(abstraction,format,mcs,channelCoding,dataLength);

perIdeal = interp1(lut(:,1),lut(:,2),binsnr,'linear','extrap');
perIdeal(perIdeal>1) = 1;

positivePerIdx = binper>0;
mse = mean((binper(positivePerIdx)'-perIdeal(positivePerIdx)).^2); % Mean Squared Error 
end