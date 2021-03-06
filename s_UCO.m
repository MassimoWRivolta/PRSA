clearvars;
clc
close all

%% UCO creation.
% UCO parameters.
ucoDuration = 60;
restDuration = 120;
baselineRR = 400;
deltaRRresponse = 400;
dbSNR = 10;
tauResponse = 5;
tauRelax = 20;
t = (0:baselineRR/1000:10*60);

[RR, u] = myUCOResponse(t, baselineRR, deltaRRresponse, ucoDuration, restDuration, tauResponse, tauRelax, dbSNR);

subplot(2, 3, [1, 2])
plot(t, u, 'LineWidth', 1.2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('UCO', 'FontSize', 12);
ylim([-0.5, 1.5]);
h = get(gcf,'CurrentAxes');
set(h,'FontSize',12);

%% PRSA and AC/DC/DR computation.

% PRSA and AC/DC/DR paramters
L = 50;
T = 10;
s = 1;
[DC, prsaDC, anchorsDC] = ACDC(RR, true, L, T, s);
[AC, prsaAC, anchorsAC] = ACDC(RR, false, L, T, s);

subplot(2, 3, [4, 5])
plot(t, RR, 'LineWidth', 1.2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('RR (ms)', 'FontSize', 12);
hold on

% Draw growing segments.
bAnchorsDC = zeros(size(RR));
bAnchorsDC(anchorsDC) = 1;
bAnchorsDC = bAnchorsDC(:);
s = find(bAnchorsDC & [1 ~bAnchorsDC(1:end-1)']');
e = find(bAnchorsDC & [~bAnchorsDC(2:end)' 1]');
intervalsDC = [s e];
RRDC = nan(size(RR));
for ii = 1:size(intervalsDC, 1)
    RRDC(intervalsDC(ii, 1):intervalsDC(ii, 2)) = RR(intervalsDC(ii, 1):intervalsDC(ii, 2));   
end
plot(t, RRDC, 'r', 'LineWidth', 1.2)

bAnchorsAC = zeros(size(RR));
bAnchorsAC(anchorsAC) = 1;
bAnchorsAC = bAnchorsAC(:);
s = find(bAnchorsAC & [1 ~bAnchorsAC(1:end-1)']');
e = find(bAnchorsAC & [~bAnchorsAC(2:end)' 1]');
intervalsAC = [s e];
RRAC = nan(size(RR));
for ii = 1:size(intervalsAC, 1)
    RRAC(intervalsAC(ii, 1):intervalsAC(ii, 2)) = RR(intervalsAC(ii, 1):intervalsAC(ii, 2));      
end
plot(t, RRAC, 'g', 'LineWidth', 1.2)

% Draw anchor points.
plot(t(anchorsDC), RR(anchorsDC), 'r.')
plot(t(anchorsAC), RR(anchorsAC), 'g.')

% Draw legend.
legend({'RR', 'DC', 'AC'}, 'FontSize', 12);

h = get(gcf,'CurrentAxes');
set(h,'FontSize',12);

subplot(2, 3, [3, 6])
plot(prsaDC, 'r', 'LineWidth', 1.2)
hold on
plot(prsaAC, 'g', 'LineWidth', 1.2)
legend({'PRSA DC', 'PRSA AC'}, 'FontSize', 12)
title(sprintf('DC: %.2f ms, AC: %.2f ms\nDC+AC: %.2f ms', DC, AC, DC + AC), 'FontSize', 12)
h = get(gcf,'CurrentAxes');
set(h,'FontSize',12);