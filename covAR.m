function rho = covAR(a, sigma)
% It estimates the autocovariance function of an autoregressive model or
% order p.
% Only the p+1 coefficients are estimated.
% 
% INPUT:
% a: row vector of the coefficient of the AR model in the form 
%     [1 a1 a2 a3 ... ap]
% sigma: standard deviation of the zero mean Gaussian noise injected into
% the linear filter.
% 
% OUTPUT:
% rho: autocovariance function.
% 
% EXAMPLE:
% rho = 0.90;
% a = [1, -2*rho*cos(pi/3), rho^2];
% sigma = 1;
% rho = covAR(a, sigma)
% 
% VERSION:
% 1.0.0 First release.
% 
% LAST UPDATE:
% 02/09/2019

P = length(a) - 1;
A1 = zeros(P + 1);
for p = 1:P + 1 
    A1(p, 1:P + 1 - p + 1) = a(p:P + 1);
end
A2 = zeros(P + 1);
for p = 2:P + 1
    A2(p, 2:p) = fliplr(a(1:p - 1));
end
A = A1 + A2;
B = zeros(P + 1, 1);
B(1, 1) = sigma^2;

rho = A\B;
end