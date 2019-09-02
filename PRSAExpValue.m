function P = PRSAExpValue(SigmaX, T, signG)
% It computes the expected value of the PRSA series of a zero-mean
% Multivariate Gaussian process defined by the covariance matrix SigmaX.
%
% INPUT:
% SigmaX: 2Lx2L covariance matrix.
% T: T value of the PRSA algorithm.
% [signG]: -1 OR +1, defining accelerations or decelerations (optional;
% default: +1).
%
% OUTPUT:
% P: 2Lx1 PRSA series.
%
% EXAMPLE 1:
% L = 50;
% a = 0.9;
% acovfunction = a.^(0:2*L-1);
% SigmaX = toeplitz(acovfunction);
% T = 10;
% P = PRSAExpValue(SigmaX, T);
% subplot(1, 2, 1)
% plot(acovfunction);
% title('Autocovariance function');
% xlabel('Lag')
% subplot(1, 2, 2)
% plot(P);
% title(sprintf('PRSA (L = %d, T = %d)', L, T))
% 
% EXAMPLE 2:
% L = 50;
% a = 0.95;
% acovfunction = (a.^(0:2*L-1)).*cos(2*pi*4*(0:2*L-1)/(2*L));
% SigmaX = toeplitz(acovfunction);
% T = 10;
% P = PRSAExpValue(SigmaX, T);
% subplot(1, 2, 1)
% plot(acovfunction);
% title('Autocovariance function');
% xlabel('Lag')
% subplot(1, 2, 2)
% plot(P);
% title(sprintf('PRSA (L = %d, T = %d)', L, T))
%
% VERSION:
% 1.0.0 First release.
%
% LAST UPDATE:
% 22/07/2019

defSignG = 1;
if(nargin >= 3)
    if(~isempty(signG) && isscalar(signG) && abs(signG) == 1)
        defSignG = signG;
    else
        warning('signG must be either +1 or -1. signG has been set to +1.');
    end
end

% L and T.
L = sqrt(numel(SigmaX))/2;
if(L ~= floor(L))
    L = floor(L);
    SigmaX = SigmaX(1:2*L, 1:2*L);
    warning('SigmaX must have an even amount of columns. L is set to floor(L) and SigmaX = SigmaX(1:2*L, 1:2*L).');
end
if(T ~= min([T, L]))
    T = min([T, L]);
    warning('T must be lesser or equal than L. It has been set to min([T, L]).');
end

% Creation of vector g for anchor point region definition.
g = zeros(2*L, 1); g(L + 1: L + T) = 1; g(L - T + 1:L) = -1;
g = defSignG*g;

% Diagonalization of SigmaX.
[U, D] = eig(SigmaX);

% Rotation with Householder transformation.
w = (g'*U*sqrt(D))';
I = eye(2*L);
eLp1 = I(:, L + 1);
v = w + norm(w)*eLp1;
if(all(v == 0))
    v = w - norm(w)*eLp1;
    H = (I - 2*(v*v')/(v'*v));
else
    H = -(I - 2*(v*v')/(v'*v));
end

% Computation of PRSA.
P = sqrt(2/pi)*U*sqrt(D)*H'*eLp1;

end