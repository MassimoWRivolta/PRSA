function [cap, prsa, anchors, stdPrsa] = ACDC(RR, deceleration, L, T, s, noAnchors)
% This function computes the acceleration and deceleration capacities using
% the PRSA approach.
% 
% INPUT:
% RR: interbeat time interval series (or other time series).
% deceleration: boolean value (if true, it will compute the deceleration).
% L: L parameter of PRSA.
% T: T parameter of PRSA.
% s: s parameter of PRSA.
% [noAnchors]: array of indexes to exclude from the anchor point's list
% (optional).
% 
% OUTPUT:
% cap: capacity value.
% prsa: PRSA vector.
% anchors: list of anchors points of RR.
% stdPrsa: standard deviation of each PRSA point.
% 
% EXAMPLE 1:
% RR = 40*randn(1, 1000)+800;
% deceleration = true;
% L = 20;
% T = 5;
% s = 3;
% [cap, prsa, anchors] = ACDC(RR, deceleration, L, T, s);
% plot(prsa);
% title(sprintf('Cap: %.3f', cap));
% ylabel('RR (ms)');
% 
% EXAMPLE 2:
% RR = 40*randn(1, 1000)+800;
% deceleration = true;
% L = 20;
% T = 5;
% s = 3;
% [cap, prsa, anchors] = ACDC(RR, deceleration, L, T, s);
% excludedAnchorPoints = 1:2:1000;
% [capNoAnchor, prsaNoAnchor, anchors] = ACDC(RR, deceleration, L, T, s, excludedAnchorPoints);
% plot(prsa);
% hold on
% plot(prsaNoAnchor);
% title(sprintf('Cap: %.3f\nCap with excluded APs: %.3f', cap, capNoAnchor));
% legend('PRSA', 'PRSA with excluded APs');
% ylabel('RR (ms)');
%
% VERSION:
% 1.0.0 First release.
%
% LAST UPDATE:
% 22/07/2019

if(nargin <= 5)
    noAnchors = [];
else
    if(~isvector(noAnchors) || ~isnumeric(noAnchors))
        noAnchors = [];
    end
end

RR = RR(:)';

% Determine anchor points.
left = filter([0, ones(1, T)/T], 1, RR);
left(1:max([T, L])) = NaN;

right = fliplr(filter(ones(1, T)/T, 1, fliplr(RR)));
right(end-max([T, L])+1:end) = NaN;

d = (right - left)*(deceleration*2 - 1);
anchors = find(d>0);
anchors = setdiff(anchors, noAnchors);

% Compute PRSA.
prsa = nan(2*L, 1);
stdPrsa = nan(2*L, 1);
for ii = 1:2*L
   prsa(ii) = nanmean(RR(anchors - L + ii - 1));
   stdPrsa(ii) = nanstd(RR(anchors - L + ii - 1));
end

% Compute average capacity.
cap = sum((prsa(L + 1 : L + min(s, L)) - prsa(L - min(s, L) + 1 : L)))/(2*min(s, L));
end