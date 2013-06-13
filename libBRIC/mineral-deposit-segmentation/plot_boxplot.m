function [p] = plot_boxplot(V, V_ref, step, Q)
[x, sidx] = sort(log10((V+V_ref)/2));
y = 2*(V-V_ref)./(V+V_ref);
y = y(sidx);
x_mat = NaN(step, ceil(length(x)/step));
x_mat(:) = x;
y_mat = NaN(step, ceil(length(x)/step));
y_mat(:) = y;
figure;
boxplot(y_mat, median(x_mat));
ylim([-2.1 2.1]);

N_Q = length(Q);
p = NaN(N_Q, 2);
for idx_Q = 1:N_Q
    xx = ones(size(x_mat, 2), 2);
    xx(:, 1) = median(x_mat);
    yy = quantile(y_mat, Q(idx_Q))';
    p(idx_Q, :) = xx\yy;
%     Tmp = mcdregres(xx(:, 1), yy, 'plots', 0, 'alpha', .8);
% 	p(idx_Q, :) = [Tmp.slope Tmp.int];
end

figure;
scatter(10.^x, y, 10, 'k');
hold on;
for idx_Q = 1:N_Q
    plot(10.^x, polyval(p(idx_Q, :), x));
end
ylim([-2.1 2.1]);

