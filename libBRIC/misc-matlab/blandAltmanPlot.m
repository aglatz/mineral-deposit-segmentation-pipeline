function [Idx, avg, SD] = blandAltmanPlot(A,B, varargin)
%reference: Y H Chan, Biostatistics 104:
%Correlational Analysis,
%Singapore Med J 2003 Vol 44(12) : 614-619
meanAB=(A+B)./2;
if ~isempty(varargin)
    difff = varargin{1};
else
    difff = (A-B)./meanAB;
end
M = isnan(difff);
meanAB(M)=[];
difff(M)=[];

% Transform data
[meanAB_sorted, meanAB_sorted_idx] = sort(meanAB);
x_old = (meanAB_sorted(:));
y_old = difff(meanAB_sorted_idx)';
x = log10(x_old);
y = y_old;

% Quantile regression
[p, stats] = quantreg(x, y, .5);
middle_sorted = polyval(p, x);
[pp] = polyfit(x, y, 1);
mean_sorted = polyval(pp, x);
[p_upper, stats] = quantreg(x, y, .95);
upper_sorted = polyval(p_upper, x);
[p_lower, stats] = quantreg(x, y, .05);
lower_sorted = polyval(p_lower, x);
if ~isempty(varargin)
    upper_sorted(upper_sorted > 1) = 1;
    lower_sorted(lower_sorted < 0) = 0;
else
    upper_sorted(upper_sorted > 2) = 2;
    lower_sorted(lower_sorted < -2) = -2;    
end
[p_upper75, stats] = quantreg(x, y, .75);
upper75_sorted = polyval(p_upper75, x);
[p_lower25, stats] = quantreg(x, y, .25);
lower25_sorted = polyval(p_lower25, x);

% Determine outliers
M = y_old >= upper_sorted | y_old <= lower_sorted;
Idx = meanAB_sorted_idx(~M);

% Plot data
scatter(meanAB(Idx), difff(Idx), 20, 'k');
hold on;
Idx = meanAB_sorted_idx(M);
scatter(meanAB(Idx), difff(Idx), 20, '+r');

% Plot regression lines
plot(meanAB_sorted, middle_sorted, 'b', 'linewidth', 2);
plot(meanAB_sorted, mean_sorted, '--k', 'linewidth', 2);
plot(meanAB_sorted, upper_sorted, ':b', 'linewidth', 1);
plot(meanAB_sorted, lower_sorted, ':b', 'linewidth', 1);
plot(meanAB_sorted, upper75_sorted, '--b', 'linewidth', 1);
plot(meanAB_sorted, lower25_sorted, '--b', 'linewidth', 1);

xlim([min(meanAB) max(meanAB)]);
xlabel('\bf Average volume in mm^3');
if ~isempty(varargin)
    ylim([-.05 1.2]);
	ylabel('\bf Jaccard index');
else
	ylim([-2.1 2.1]);
	ylabel('\bf (Vol-Vol_{ref}) per average volume');
end
h = findobj(gcf, 'type', 'line');
legend([h(end) h(2) h(4)], 'Median','50% Range', '90% Range');