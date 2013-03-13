function [] = plot_fitprobdistunivparamhist(Hist, Dist)
% Plots the histogram returned by fitprobdistunivparam().
% INPUTS: Hist - the histogram
%         Dist - the distributions that were fitted to the data
bar(Hist);
Lim = axis;
Lim(end) = 100;
axis(Lim);
set(gca, 'XTickLabel', [{'N/A'}, Dist]);
xlabel('\bf Distribution type');
set(gca, 'YTick', 0:20:100);
set(gca, 'YTickLabel', {'0%', '20%', '40%', '60%', '80%', '100%'});
ylabel('\bf Occurence in percent');
