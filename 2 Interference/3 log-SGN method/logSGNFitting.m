function xBest = logSGNFitting(gammaEffLinear)
logGammaEffLinearMix = log(gammaEffLinear);
n = length(logGammaEffLinearMix);
M = mean(logGammaEffLinearMix);
S = std(logGammaEffLinearMix);
xBest = zeros(1,4);
fBest = -10e10;
randRange1 = 1;
randRange2 = 0.5;
x0 = [M,S, -randRange1 + (2*randRange1)*rand, (2*randRange2)*rand];
A = [0, 0, 0, -1; 0, -1, 0, 0];
b = [0, 0];
for iter = 1: 100
    LikelihoodFun = @(x) -(-n/2* log(pi * x(2)^2/2) - 1/2/x(2)^2*sum((logGammaEffLinearMix-x(1)).^2) + sum(log(normcdf(x(3).*(logGammaEffLinearMix-x(1))./sqrt(x(2)^2+ x(4).*(logGammaEffLinearMix-x(1)).^2)))));
    [x, fval] = fmincon(LikelihoodFun,x0,A,b);
    if fval>fBest
        fBest = fval;
        xBest = x;
        x0 = xBest;
    end
end
end