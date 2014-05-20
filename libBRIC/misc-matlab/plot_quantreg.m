function [idx_oli, p_med] = plot_quantreg(x_in, y_in, y_min, y_max)
% Transform data
[x, idx_sorted] = sort(x_in(:));
y = y_in(idx_sorted); y = y(:);

% Linear model
p_mean = robustfit(x, y);
y_mean = polyval([p_mean(2) p_mean(1)], x);

% Quantile regression
Q = [.05 .25 .5 .75 .95];
N_Q = length(Q);
y_quant = zeros(length(x), N_Q);
p_med = [NaN; NaN];
for idx_q = 1:N_Q
    p = quantreg(x, y, Q(idx_q));
    if Q(idx_q) == .5
        p_med = p;
    end
    y_quant(:, idx_q) = polyval(p, x);
end

% Determine outliers
M = y >= y_quant(:, end) | y <= y_quant(:, 1);
idx_oli = idx_sorted(M);

% limits
for idx_q = 1:N_Q
    y_quant(y_quant(:, idx_q) > y_max, idx_q) = NaN;
    y_quant(y_quant(:, idx_q) < y_min, idx_q) = NaN;
end

% Plot data
scatter(x(~M), y(~M), 20, 'k', 'filled');
hold on;
scatter(x(M), y(M), 20, '+r');

% Plot mean
% plot(x, y_mean, '--k', 'linewidth', 1);

% Plot regression lines
F = {':b', '--b', 'b', '--b', ':b'};
for idx_q = 1:N_Q
    plot(x, y_quant(:, idx_q), F{idx_q}, 'linewidth', 2);
end
