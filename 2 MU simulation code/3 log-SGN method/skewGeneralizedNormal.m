function X = skewGeneralizedNormal(epsilo, omega, lamdbda1, lambda2)
alpha = sqrt(lambda2) * normrnd(0,1) + lamdbda1;
X = skewNormal(epsilo, omega, alpha);
end