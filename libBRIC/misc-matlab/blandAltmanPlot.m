function [idx_oli] = blandAltmanPlot(A,B, varargin)
%reference: Y H Chan, Biostatistics 104:
%Correlational Analysis,
%Singapore Med J 2003 Vol 44(12) : 614-619
meanAB=(A+B)./2;
if ~isempty(varargin)
    difff = varargin{1};
else
    difff = (A-B)./meanAB;
end

% Remove NaNs
M = isnan(difff);
meanAB(M)=[];
difff(M)=[];

meanAB_log10 = log10(meanAB);

% Plot scatter plut quantile regression lines
[idx_oli] = plot_quantreg(meanAB_log10, difff);

% Annotate axis
xtick = unique(round(get(gca, 'XTick')));
xlim([xtick(1)-0.05 xtick(end)+0.05]);
set(gca, 'XTick', xtick);
set(gca, 'XTickLabel', num2str(10.^(xtick')));

xlabel('\bf Average volume in mm^3');
if ~isempty(varargin)
    ylim([-.05 1.2]);
	ylabel('\bf Jaccard index');
else
	ylim([-2.1 2.1]);
	ylabel('\bf (Vol-Vol_{ref}) per average volume');
end
% h = findobj(gcf, 'type', 'line');
% legend([h(end) h(2) h(4)], 'Median','50% Range', '90% Range');