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
plot(t(anchorsDC), RR(anchorsDC), 'r*');
plot(t(anchorsAC), RR(anchorsAC), 'g*');
legend({'RR', 'Anchors DC', 'Anchors AC'}, 'FontSize', 12)
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