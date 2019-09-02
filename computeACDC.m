function cap = computeACDC(prsa, s)
% It computes the average capacity of the input PRSA series.
% 
% INPUT:
% prsa: vector of the PRSA series.
% s: s parameter.
%
% OUTPUT:
% cap: average capacity.
%
% VERSION:
% 1.0.0 First release.
%
% LAST UPDATE:
% 02/09/2019

if(mod(length(prsa), 2) ~= 0)
   error('PRSA length must be even.');   
end

L = length(prsa)/2;

% Average capacity computation.
cap = sum((prsa(L + 1 : L + min(s, L)) - prsa(L - min(s, L) + 1 : L)))/(2*min(s, L));

end