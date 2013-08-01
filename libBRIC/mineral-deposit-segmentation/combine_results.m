close all; clear all; clc

Name = 'subjects_98';
Pfx = 'ad_100713_lnorm0.6';
N_cpus = 6;
AdaptiveFlag = true;
RoiLabelTable = {[13 11 12 14]};
% [Ret, Subjects, Over] = segment_us_mp([Name '.xls'], ...
%                                       RoiLabelTable, N_cpus, ...
%                                       'ThreshFactor', [1 0], ...
%                                       'AdaptiveFlag', AdaptiveFlag, ...
%                                       'SaveMaskFlag', true, ...
%                                       'IntvarP', 0.4);
% save([Name '_' Pfx '.mat']);

[J, D, V, V_ref, Subjects, FPC, FPV, FPSubjects, Sens, Spec, Conf] = load_matdata([Name '_' Pfx '.mat'], 0);

% AB = reshape([out.P], 2, 10);
% fprintf('-- alpha --');
% quantile(AB(1, :), [.1 .5 .9])
% % bootci(10000, {@median, AB(1, :)}, 'alpha', 0.05, 'type', 'bca')
% fprintf('-- beta --');
% quantile(AB(2, :), [.1 .5 .9])
% % bootci(10000, {@median, AB(2, :)}, 'alpha', 0.05, 'type', 'bca')

fprintf('-- Jaccard --');
Q = [.25 .5 .75];
quantile(J, Q)
% bootci(10000, {@median, J(M_vol_ref)}, 'alpha', 0.05, 'type', 'bca')
fprintf('-- Dice --');
quantile(D, Q)
% bootci(10000, {@median, D(M_vol_ref)}, 'alpha', 0.05, 'type', 'bca')
fprintf('-- Volume --');
quantile(V, Q)
% bootci(10000, {@median, V(M_vol_ref)}, 'alpha', 0.05, 'type', 'bca')
fprintf('-- FP count --');
FPC
fprintf('-- FP volume --');
quantile(FPV, Q)
% bootci(10000, {@median, V(M_vol & ~M_vol_ref)}, 'alpha', 0.05, 'type', 'bca')
fprintf('-- Sens/Spec --');
[Sens Spec]
fprintf('-- Confusion --');
Conf
% fprintf('-- Param --');
% quantile(P, Q)

save_xls('FP', FPSubjects);

figure;
Idx = blandAltmanPlot(V, V_ref)';
set(gca, 'XScale', 'log');
save_xls('BA_OLI', Subjects(Idx));
set(gcf, 'color', 'white')

% Testing
% plot_boxplot(V, V_ref, 15, [.05, .25, .5, .75, .95]);
% set(gca, 'XScale', 'log');
% set(gcf, 'color', 'white');

figure;
Idx = blandAltmanPlot(V, V_ref, J)';
set(gca, 'XScale', 'log');
save_xls('MBA_OLI', Subjects(Idx));
set(gcf, 'color', 'white')

save('forR', 'V', 'V_ref', 'J', '-v6');
% system('R CMD BATCH myqr.R');

% M_over = sum(isnan(Over),2)==0;
% [R,P] = corrcoef(Over(M_over, :))
