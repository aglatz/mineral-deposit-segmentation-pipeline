function [H] = summary_fitprobdistunivparam(Bfit, Dist_label, ti, F, ICV)
% Plots the average distribution shape of each fitted distribution.
% INPUTS: Bfit - Cell array with the results from fitprobdistunivparam()
%         Dist_label - Cell array with distribution names that were
%                      fitted to the data
%         ti - title of the plot
%         F - voxel size
%         ICV - Intracranial volume of every subject in Bfit
% RETURNS: H - figure handle
%
N_lab = length(Dist_label);
H = figure;
for idx_lab = 1:N_lab
    M = [Bfit{:, 1}] == idx_lab;
    if sum(M) > 0
        Tmp = Bfit(M, 2);
        %Tmp1 = I(M);
        Params = zeros(sum(M), length(Tmp{1}.Params));
        for idx = 1:sum(M)
            Params(idx, :) = Tmp{idx}.Params;
            
            %plot_mydata(Tmp1{idx}, I_name, Dist_label, Dist_label{idx_lab});
        end
        fprintf('- %s ------------------------------------------\n', Dist_label{idx_lab});
        N = sum(M)
        Q = quantile(Params, [.25 .5 .75], 1)
        V = quantile([Bfit{M, 3}]./ICV(M)*1e6, [.25 .5 .75])
    else
        N = 0;
        Q = [];
        V = [0 0 0];
    end
	subplot(1, N_lab, idx_lab); hold on;
    N_Q = size(Q, 1);
    X = 0.01:0.01:0.99;
    %Map = gray(4);
    for idx_q = 1:N_Q
        PD = ProbDistUnivParam(Dist_label{idx_lab}, Q(idx_q, :));
        Y = PD.pdf(X);
        if idx_q == 1 || idx_q == 3
            plot(X, Y, '--k', 'LineWidth', 1);
        else
            plot(X, Y, 'k', 'LineWidth', 2);
        end
    end
    xlabel(['\bf Transformed ' ti ' intensities']);
    ylabel('\bf Probability density');
    Lim = [0 1 0 3.1];
    axis(Lim);
    if size(Q, 2) > 0
        Param_str = sprintf('(p%d=%0.1f', 1, Q(2, 1));
        for idx_p = 2:size(Q, 2)
            Param_str = [Param_str sprintf(', p%d=%0.1f', idx_p, Q(2, idx_p))];
        end
        Param_str = [Param_str ')'];
    else
        Param_str = '';
    end
    text(Lim(2)*0.1, Lim(end)*0.9, sprintf(['\\bf %s %s\n N=%d\n ' ... 
         'V^{norm}=%0.1fppm (%0.1f...%0.1f)ppm'], ...
         Dist_label{idx_lab}, Param_str, N, ...
         get_volume(V(2), F), get_volume(V(1), F), get_volume(V(3), F)));
end
%legend('25^{th} percentile', '50^{th} percentile', '75^{th} percentile', 'location', 'Best');
fprintf('===============================================\n\n');


function [] = plot_mydata(I, I_name, Dist_label, Dist_best)
X_dbl = double(I);
P = polyfit([min(X_dbl) max(X_dbl)], [.99 .01], 1);
X_tf = polyval(P, X_dbl);
X = 0.01:0.01:0.99;
N_dist = length(Dist_label);
Col = hsv(6);
Max_val = zeros(N_dist, 2);
Y = cell(N_dist, 1);
H1 = figure;
for idx = 1:N_dist
    PD = fitdist(X_tf(:), Dist_label{idx});
    Y1 = PD.pdf(X);
    idx_max = find(Y1==max(Y1));
    Max_val(idx, :) = [X(idx_max(1)) Y1(idx_max(1))];
    subplot(3,2,idx*2);
    hold on; qqplot(X_tf(:), PD);
    Y{idx} = Y1;
end
Max_val_med = median(Max_val(1:2, :), 1);
subplot(3, 2, [1 3]);
plot_hist(X_tf, Max_val_med(2), Max_val_med(1), 'b', 1);
for idx = 1:N_dist
    hold on; plot(X, Y{idx}, 'color', Col(idx*2, :), 'LineWidth', 1);
end
Lim = axis;
Lim(end) = Lim(end)+.4;
axis(Lim);
legend([{'histogram', 'kerneldensity'} Dist_label]);
title(Dist_best);
save_ps_figure(I_name, H1);
