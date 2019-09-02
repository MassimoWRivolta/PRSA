function [out, sigma] = acovfun(a, L)
% It builds the autocovariance function of length L and its corresponding
% covariance matrix of an autoregressive process defined by the vector
% parameter a.
% 
% It assumes the variance of the input noise to be 1.
% 
% It is worth noting that the maximum lag considered is L-1.
% 
% INPUT:
% a: coefficients of the autoregressive process.
% L: length of the autocovariance function and dimension of the covariance
% matrix.
% 
% OUTPUT:
% out: array with the autocovariance function of length L.
% sigma: LxL covariance matrix.
% 
% EXAMPLE:
% rho = 0.90;
% a = [1, -2*rho*cos(pi/3), rho^2];
% L = 50;
% [out, sigma] = acovfun(a, L);
% plot(0:L-1, out);
% xlabel('Lag')
% ylabel('Amplitude')
% title('Autocovariance function')
% 
% DEPENDENCIES:
% covAR.m
% 
% VERSION:
% 1.0.0 First release.
% 
% LAST UPDATE:
% 02/09/2019

if(~exist('covAR.m', 'file'))
    error('acovfun requires the function covAR.m');
end

rho = covAR(a, 1);
rho = rho(:)';

% Compute autocovariance function up to lag L-1.
N = length(rho);
out = zeros(1, L);
out(1:min([N, L])) = rho(1:min([N, L]));
for nn = N + 1:L
   out(nn) = -sum(fliplr(a(2:end)).*out(nn-N+1:nn-1));
end

sigma = zeros(L, L);

% Build covariance matrix.
for rr = 1:L
    for cc = rr:L
        sigma(rr, cc) = out(abs(rr - cc) + 1);
        sigma(cc, rr) = sigma(rr, cc);
    end
end

end