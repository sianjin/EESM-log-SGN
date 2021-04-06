% Generate SGN random variable 
% shown in Algorithm. 2 in the IEEE TCOM paper:
% "Efficient PHY Layer Abstraction for Fast Simulations in Complex 
% System Environments"
function X = skewGeneralizedNormal(epsilo, omega, lamdbda1, lambda2)
alpha = sqrt(lambda2) * normrnd(0,1) + lamdbda1;
X = skewNormal(epsilo, omega, alpha);
end