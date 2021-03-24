function plotSINR(sinr,Hsoi,Psoi,Hint,Pint,N0,k)
% plotRBIR Example helper function

%   Copyright 2019 The MathWorks, Inc.

% Reshape so all combinations of tx and rx antenna are in the 2nd dimension
Hsoi = reshape(Hsoi,size(Hsoi,1),[],1,1);
Hint = reshape(Hint,size(Hint,1),[],1,1);

figure;
subplot(211);
% Plot signal of interest
hs = plot(k,10*log10(abs(Hsoi).^2)+Psoi);
for i = 2:numel(hs) % Color of all tx-rx antenna combination is the same
    hs(i).Color = hs(1).Color;
end
hold on;
% Plot interfering interest
hi = plot(k,10*log10(abs(Hint).^2)+Pint);
for i = 2:numel(hi) % Color of all tx-rx antenna combination is the same
    hi(i).Color = hi(1).Color;
end

% Plot noise floor
hn = plot(k,repmat(N0,size(k)));
xlabel('Subcarrier');
ylabel('dBm');
title('Per-subcarrier receive power');
grid on;
legend([hs(1) hi(1) hn],'Signal of interest','Interferer','Noise floor','location','best');

subplot(212);
% Plot calculated SINR per spatial stream
plot(k,sinr); 
xlabel('Subcarrier');
ylabel('dB');
title('Post-equalizer SINR');
grid on;

end