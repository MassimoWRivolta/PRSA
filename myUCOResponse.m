function [RR, u] = myUCOResponse(t, baselineRR, deltaResponseRR, ucoDuration, restDuration, tauResponse, tauRelax, dbSNR)
% myUCOResponse generates the response of an exponential response to a
% square wave.
%
% INPUT:
% t: time vector (s).
% baselineRR: average baseline RR (ms).
% deltaResponseRR: maximum change in RR (ms).
% ucoDuration: duration (s).
% restDuration: rest duration, between two UCOs (s).
% tauResponse: time constant to the UCO response.
% tauRelax: time constant when UCO finishes.
% dbSNR: Signal-to-Noise ratio in db.
%
% OUTPUT:
% RR: interbeat interval time series (ms).
% u: UCO signal.
%
% EXAMPLE:
% ucoDuration = 60; 
% restDuration = 60;
% baselineRR = 400;
% deltaRRresponse = 400;
% dbSNR = 10;
% tauResponse = 5;
% tauRelax = 20;
% t = (0:baselineRR/1000:10*60);
% [RR, u] = myUCOResponse(t, baselineRR, deltaRRresponse, ucoDuration, restDuration, tauResponse, tauRelax, dbSNR);
% subplot(2, 1, 1)
% plot(t, u);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('UCO');
% ylim([-0.5, 1.5]);
% subplot(2, 1, 2)
% plot(t, RR);
% xlabel('Time (s)');
% ylabel('RR (ms)');
%
% VERSION:
% 1.0.0 First release.
%
% LAST UPDATE:
% 02/09/2019

param(1) = tauResponse;
param(2) = tauRelax;
param(3) = ucoDuration;
param(4) = restDuration;

x0 = 0;
solution  =  ode23(@(tt, x) myUCODerivative(tt, x, param), [t(1), t(end)], x0);
data = deval(solution, t);
max_Y = max(data);
min_Y = min(data);
data = (data - min_Y)/(max_Y - min_Y)*deltaResponseRR + baselineRR;
stdNoise = sqrt(var(data)/10^(dbSNR/10));

RR = data + stdNoise*randn(size(t));

if(nargout == 2)
   u = myUCO(t, ucoDuration, restDuration); 
end
end

function dydt = myUCODerivative(t, y, param)
% myUCODerivative defines the ODE for an exponential response to a square
% wave.
%
% INPUT:
% t: time vector (unit ut).
% y: variable value (unit uy).
% param: 1x4 vector of parameters
%   param(1) = tauResponse (ut)
%   param(2) = tauRelax (ut)
%   param(3) = ucoDuration (ut)
%   param(4) = restDuration (ut)
%
% OUTPUT:
% dydt: derivative of y with respect to t (uy/ut).
%
% EXAMPLE:
% t = 0:1/1000:10*60; % (s).
% param(1) = 5;
% param(2) = 10;
% param(3) = 60; % (s).
% param(4) =  30; % (s).
% y = randn(size(t));
% dydt = myUCOResponse(t, y, param);
% plot(t, dydt);
% xlabel('Time (s)');
% ylabel('dy/dt');
%
% VERSION:
% 1.0.0 First release.
%
% LAST UPDATE:
% 02/09/2019

tau1 = param(1);
tau2 = param(2);
ucoDuration = param(3);
restDuration = param(4);

u = myUCO(t, ucoDuration, restDuration);

dydt = -(1/tau2)*y - (1/tau1 - 1/tau2)*(u.*y) + u/tau1;
end