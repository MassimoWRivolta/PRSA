function P = PRSAExpValueAR(a, sigma, L, T, signG)
% It computes the expected value of the PRSA series of an autoregressive
% model defined by the coefficient vector [1 a1 a2 a3 ... ap].
% 
% INPUT:
% a: vector of the coefficient of the AR model.
% sigma: standard deviation of the input white Guassian noise.
% L: L parameter of the PRSA algorithm.
% T: T parameter of the PRSA algorithm.
% [signG]: -1 OR +1, defining accelerations or decelerations (optional;
% default: +1).
%
% OUTPUT:
% P: 2Lx1 PRSA series.
%
% EXAMPLE:
% rho = 0.90;
% a = [1, -2*rho*cos(pi/3), rho^2];
% sigma = 1;
% L = 50;
% T = 10;
% P = PRSAExpValueAR(a, sigma, L, T);
% plot(P);
% title(sprintf('PRSA (L = %d, T = %d)', L, T))
%
% DEPENDENCIES:
% acovfun.m
% 
% VERSION:
% 1.0.0 First release.
%
% LAST UPDATE:
% 02/09/2019

if(~exist('acovfun.m', 'file'))
    error('PRSAExpValueAR requires the function acovfun.m');
end

defSignG = 1;
if(nargin >= 5)
    if(~isempty(signG) && isscalar(signG) && abs(signG) == 1)
        defSignG = signG;
    else
        warning('signG must be either +1 or -1. signG has been set to +1.');
    end
end

if(L ~= max([2, abs(L)]))
    L = max([2, abs(L)]);
    warning('L must be greater or equal than 2. It has been set to 2.');
end

[~, SigmaX] = acovfun(a, 2*L);
SigmaX = (sigma^2)*SigmaX;

% L and T.
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