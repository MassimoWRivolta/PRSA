clearvars;
close all
clc

A = 0.99;
N = 500;
w = linspace(0, pi, N);

L = 50;

T = [1, 1, 1, 1, 2, 10, 2, 10, 10, 2];
s = [1, 2, 10, 5, 2, 10, 1, 1, 2, 10];
nPairs = length(T);

cap = zeros(N, nPairs);
for pp = 1:nPairs
    for ii = 1:N
        a = [1, -2*A*cos(w(ii)), A^2];
        
        [~, SigmaX] = acovfun(a, 2*L);
        var0 = SigmaX(1, 1);
        SigmaX = SigmaX/var0;
        
        P = PRSAExpValue(SigmaX, T(pp));
        
        cap(ii, pp) = computeACDC(P, s(pp));
    end
end

colors = [mod(1:nPairs, 3)' == 2, 0.7*circshift(mod(1:nPairs, 3)' == 2, 1), circshift(mod(1:nPairs, 3)' == 2, 2)];

%% T = 1, s = ...

iToConsider = [find(T == 1 & s == 1), find(T == 1 & s == 2), find(T == 1 & s == 5), find(T == 1 & s == 10)];

figure
hold on
legStr = {};

% Plot frequency response.
for pp = 1:length(iToConsider)
    plot(w/pi, cap(:, iToConsider(pp)), '.-', 'Color', colors(pp, :), 'LineWidth', 1.2);
    legStr = cat(2, legStr, sprintf('s = %d', s(iToConsider(pp))));
end
xlabel('Normalized frequency', 'FontSize', 12);
ylabel('Average capacity (a.u.)', 'FontSize', 12);
title('T = 1', 'FontSize', 12);
legend(legStr, 'FontSize', 12, 'Location', 'best');

% Plot the maximum theoretical value.
for pp = 1:length(iToConsider)
    [~, iPred] = min(abs(w - 0.371*2*pi/s(iToConsider(pp))));
    plot(w(iPred)/pi, cap(iPred, iToConsider(pp)), 'X', 'MarkerSize', 12, 'Color', colors(pp, :), 'LineWidth', 1.2)
end

% Plot the zeros of (z^s - 1)^2
for pp = 1:length(iToConsider)
    c = zeros(1, s(iToConsider(pp)) + 1);
    c(1) = 1;
    c(end) = -1;
    r = angle(roots(conv(c, c)))/pi;    
    r = sort(r(1:2:end));
    
    plot(r, 0*(pp - 1)*ones(size(r)), 's', 'MarkerSize', 8, 'Color', colors(pp, :), 'MarkerFaceColor', colors(pp, :));
end

h = get(gcf,'CurrentAxes');
set(h,'FontSize',12);

%% T = s, s = ...

iToConsider = [find(T == 1 & s == 1), find(T == 2 & s == 2), find(T == 10 & s == 10)];

figure
hold on
legStr = {};

% Plot frequency response.
for pp = 1:length(iToConsider)
    plot(w/pi, cap(:, iToConsider(pp)), '.-', 'Color', colors(pp, :), 'LineWidth', 1.2);
    legStr = cat(2, legStr, sprintf('s = %d', s(iToConsider(pp))));
end
xlabel('Normalized frequency', 'FontSize', 12);
ylabel('Average capacity (a.u.)', 'FontSize', 12);
title('T = s', 'FontSize', 12);
legend(legStr, 'FontSize', 12, 'Location', 'best');

% Plot the maximum theoretical value.
for pp = 1:length(iToConsider)
    [~, iPred] = min(abs(w - 0.371*2*pi/s(iToConsider(pp))));
    plot(w(iPred)/pi, cap(iPred, iToConsider(pp)), 'X', 'MarkerSize', 12, 'Color', colors(pp, :), 'LineWidth', 1.2)
end

% Plot the zeros of (z^s - 1)^2
for pp = 1:length(iToConsider)
    c = zeros(1, s(iToConsider(pp)) + 1);
    c(1) = 1;
    c(end) = -1;
    r = angle(roots(conv(c, c)))/pi;    
    r = sort(r(1:2:end));
    
    plot(r, 0*(pp - 1)*ones(size(r)), 's', 'MarkerSize', 8, 'Color', colors(pp, :), 'MarkerFaceColor', colors(pp, :));
end

h = get(gcf,'CurrentAxes');
set(h,'FontSize',12);

%% T = ..., s = 1

iToConsider = [find(T == 1 & s == 1), find(T == 2 & s == 1), find(T == 10 & s == 1)];

figure
hold on
legStr = {};

% Plot frequency response.
for pp = 1:length(iToConsider)
    plot(w/pi, cap(:, iToConsider(pp)), '.-', 'Color', colors(pp, :), 'LineWidth', 1.2);
    legStr = cat(2, legStr, sprintf('T = %d', T(iToConsider(pp))));
end
xlabel('Normalized frequency', 'FontSize', 12);
ylabel('Average capacity (a.u.)', 'FontSize', 12);
title('s = 1', 'FontSize', 12);
legend(legStr, 'FontSize', 12, 'Location', 'best');

% Plot the maximum theoretical value.
for pp = 1:length(iToConsider)
    [~, iPred] = min(abs(w - 0.371*2*pi/s(iToConsider(pp))));
    plot(w(iPred)/pi, cap(iPred, iToConsider(pp)), 'X', 'MarkerSize', 12, 'Color', colors(pp, :), 'LineWidth', 1.2)
end

% Plot the zeros of (z^T - 1)^2
for pp = 1:length(iToConsider)
    c = zeros(1, T(iToConsider(pp)) + 1);
    c(1) = 1;
    c(end) = -1;
    r = angle(roots(conv(c, c)))/pi;    
    r = sort(r(1:2:end));
    
    plot(r, 0*(pp - 1)*ones(size(r)), 'O', 'MarkerSize', 8, 'Color', colors(pp, :), 'MarkerFaceColor', colors(pp, :));
end

h = get(gcf,'CurrentAxes');
set(h,'FontSize',12);

%% T = ..., s = 2

iToConsider = [find(T == 1 & s == 2), find(T == 2 & s == 2), find(T == 10 & s == 2)];

figure
hold on
legStr = {};

% Plot frequency response.
for pp = 1:length(iToConsider)
    plot(w/pi, cap(:, iToConsider(pp)), '.-', 'Color', colors(pp, :), 'LineWidth', 1.2);
    legStr = cat(2, legStr, sprintf('T = %d', T(iToConsider(pp))));
end
xlabel('Normalized frequency', 'FontSize', 12);
ylabel('Average capacity (a.u.)', 'FontSize', 12);
title('s = 2', 'FontSize', 12);
legend(legStr, 'FontSize', 12, 'Location', 'best');

% Plot the maximum theoretical value.
for pp = 1:length(iToConsider)
    [~, iPred] = min(abs(w - 0.371*2*pi/s(iToConsider(pp))));
    plot(w(iPred)/pi, cap(iPred, iToConsider(pp)), 'X', 'MarkerSize', 12, 'Color', colors(pp, :), 'LineWidth', 1.2)
end

% Plot the zeros of (z^T - 1)^2
for pp = 1:length(iToConsider)
    c = zeros(1, T(iToConsider(pp)) + 1);
    c(1) = 1;
    c(end) = -1;
    r = angle(roots(conv(c, c)))/pi;    
    r = sort(r(1:2:end));
    
    plot(r, 0*(pp - 1)*ones(size(r)), 'O', 'MarkerSize', 8, 'Color', colors(pp, :), 'MarkerFaceColor', colors(pp, :));
end

%% T = ..., s = 10

iToConsider = [find(T == 1 & s == 10), find(T == 2 & s == 10), find(T == 10 & s == 10)];

figure
hold on
legStr = {};

% Plot frequency response.
for pp = 1:length(iToConsider)
    plot(w/pi, cap(:, iToConsider(pp)), '.-', 'Color', colors(pp, :), 'LineWidth', 1.2);
    legStr = cat(2, legStr, sprintf('T = %d', T(iToConsider(pp))));
end
xlabel('Normalized frequency', 'FontSize', 12);
ylabel('Average capacity (a.u.)', 'FontSize', 12);
title('s = 10', 'FontSize', 12);
legend(legStr, 'FontSize', 12, 'Location', 'best');

% Plot the maximum theoretical value.
for pp = 1:length(iToConsider)
    [~, iPred] = min(abs(w - 0.371*2*pi/s(iToConsider(pp))));
    plot(w(iPred)/pi, cap(iPred, iToConsider(pp)), 'X', 'MarkerSize', 12, 'Color', colors(pp, :), 'LineWidth', 1.2)
end

% Plot the zeros of (z^T - 1)^2
for pp = 1:length(iToConsider)
    c = zeros(1, T(iToConsider(pp)) + 1);
    c(1) = 1;
    c(end) = -1;
    r = angle(roots(conv(c, c)))/pi;    
    r = sort(r(1:2:end));
    
    plot(r, 0*(pp - 1)*ones(size(r)), 'O', 'MarkerSize', 8, 'Color', colors(pp, :), 'MarkerFaceColor', colors(pp, :));
end

h = get(gcf,'CurrentAxes');
set(h,'FontSize',12);