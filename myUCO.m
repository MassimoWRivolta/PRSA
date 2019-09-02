function u = myUCO(t, ucoDuration, restDuration)
% myUCO generates a periodic square wave.
% The first UCO is set at t=0.
% 
% INPUT:
% t: time vector (s).
% ucoDuration: duration (s).
% restDuration: rest duration, between two UCOs (s).
% 
% OUTPUT:
% u: size(t) vector (it is 1 when the UCO is present).
% 
% EXAMPLE:
% t = 0:1/1000:10*60; % (s).
% ucoDuration = 60; % (s).
% restDuration = 30; % (s).
% u = myUCO(t, ucoDuration, restDuration);
% plot(t, u);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('UCO');
% ylim([-0.5, 1.5]);
%
% VERSION:
% 1.0.0 First release.
%
% LAST UPDATE:
% 02/09/2019

u = zeros(size(t));
u((t >= 0) & (mod(t, ucoDuration + restDuration) < ucoDuration)) = 1;

end

