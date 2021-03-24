function mse = sinrModelFittingMse(effSnrSigSGNLinear,effSnrIntSGNLinear,sinrEffdB,theta)
% Discretize true effective SINR into multiple bins with step size 0.25dB
finiteSinrEffIdx = isfinite(sinrEffdB);
sinrEffdB = sinrEffdB(finiteSinrEffIdx);
edges = floor(min(sinrEffdB(:))-20):0.25:ceil(max(sinrEffdB(:))+20);
binsinrDiscrete = discretize(sinrEffdB(:),edges);
pmfSinr = zeros(numel(edges)-1,1);
for i = 1:numel(edges)-1
    pmfSinr(i) = sum(binsinrDiscrete==i)/length(sinrEffdB);
end

% Obtain accumulative effective INR
effSnrIntSGNLinearAccum = effSnrIntSGNLinear{1};
for iter = 2:numel(effSnrIntSGNLinear)
    effSnrIntSGNLinearAccum = effSnrIntSGNLinearAccum + effSnrIntSGNLinear{iter};
end
% Obtain modeled effective SINR
sinrEffModelLinear = effSnrSigSGNLinear./(1+theta*effSnrIntSGNLinearAccum);
sinrEffModeldB = 10*log10(sinrEffModelLinear);
finiteSinrEffModelIdx = isfinite(sinrEffModeldB);
sinrEffModeldB = sinrEffModeldB(finiteSinrEffModelIdx);
binsinrModelDiscrete = discretize(sinrEffModeldB(:),edges);
pmfSinrModel = zeros(numel(edges)-1,1);
for i = 1:numel(edges)-1
    pmfSinrModel(i) = sum(binsinrModelDiscrete==i)/length(sinrEffModeldB);
end

mse = mean((pmfSinr-pmfSinrModel).^2); % Mean Squared Error 
end