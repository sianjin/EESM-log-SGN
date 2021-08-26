% Optimize log-SGN parameters by fitting the effective SINR historgram
% using maximum likelihood estimation
% Output: xBest is the vector of the 4 optimized log-SGN parameters
function xBest = logSGNFitting(gammaEffLinear)
logGammaEffLinear = log(gammaEffLinear);
% Initialize log-SGN parameters
n = length(logGammaEffLinear);
M = mean(logGammaEffLinear);
S = std(logGammaEffLinear);
xBest = zeros(1,4);
fBest = -10e10;
randRange1 = 2;
randRange2 = 2;
x0 = [M,S, -randRange1 + (2*randRange1)*rand, (2*randRange2)*rand]; % Initial log-SGN parameters
A = [0, 0, 0, -1; 0, -1, 0, 0];
b = [0, 0];
% Optimize log-SGN parameters using maximum likelihood estimation
for iter = 1: 100
    LikelihoodFun = @(x) -(-n/2* log(pi * x(2)^2/2) - 1/2/x(2)^2*sum((logGammaEffLinear-x(1)).^2) + sum(log(normcdf(x(3).*(logGammaEffLinear-x(1))./sqrt(x(2)^2+ x(4).*(logGammaEffLinear-x(1)).^2)))));
    [x, fval] = fmincon(LikelihoodFun,x0,A,b);
    if fval>fBest
        fBest = fval;
        xBest = x;
        x0 = xBest;
    end
end
end