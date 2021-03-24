function X = skewNormal(epsilo, omega, alpha)
m = sqrt((1 + alpha^2)/2) * epsilo;
U1 = normrnd(m,omega);
U2 = normrnd(m,omega);
U = max (U1, U2);
V = min (U1, U2);
lamda1 = (1+alpha)/sqrt(2*(1+alpha^2));
lamda2 = (1-alpha)/sqrt(2*(1+alpha^2));
X = lamda1 * U + lamda2 * V;
end