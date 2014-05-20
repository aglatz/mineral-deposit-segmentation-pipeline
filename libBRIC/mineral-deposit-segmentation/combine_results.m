function [Avg, SD, Conf] = combine_results()

Name = 'subjects_98';
Pfx = 'ad_t1w_mineral_test';
AdaptiveFlag = true;
RoiLabelTable = {[13, 11, 12, 14]};
N_cpus = 6;
% [xopt] = segment_us_refine([Name '.xls'], ...
%                            RoiLabelTable, ...
%                            N_cpus, ...
%                            AdaptiveFlag);
[Ret, Subjects, Over] = segment_us_mp([Name '.xls'], ...
                                      RoiLabelTable, N_cpus, ...
                                      'ThreshFactor', [1 0], ...
                                      'AdaptiveFlag', AdaptiveFlag, ...
                                      'SaveMaskFlag', true, ...
                                      'ReportName', 'class', ...
                                      'N_gre', 1, ...
                                      'CNR_thr', 0, ...
                                      'phypo_thr', 0.1, ...
                                      'intvar_thr', 0.8);
save([Name '_' Pfx '.mat']);

[Subjects, J, D, V, V_ref] = load_matdata([Name '_' Pfx '.mat'], 0);

fprintf('-- Jaccard (All) --');
i=1; [Avg(i), SD(i)] = get_avg_sd(J)

fprintf('-- Dice (All) --');
i=2; [Avg(i), SD(i)] = get_avg_sd(D)

fprintf('-- Volume (All) --');
i=3; [Avg(i), SD(i), Q] = get_avg_sd(V)

fprintf('-- Reference volume (All) --');
i=4; [Avg(i), SD(i), Q] = get_avg_sd(V_ref)

% fprintf('-- IntVarP (All) --');
% quantile(IntVarP, Q)

M_vol = V > 0;
M_vol_ref = V_ref > 0;
M = M_vol | M_vol_ref; % & J > 0;

fprintf('-- Jaccard --');
i=5; [Avg(i), SD(i)] = get_avg_sd(J(M))

fprintf('-- Dice --');
i=6; [Avg(i), SD(i)] = get_avg_sd(D(M))

fprintf('-- Volume --');
i=7; [Avg(i), SD(i), Q] = get_avg_sd(V(M))

fprintf('-- Reference volume --');
i=8; [Avg(i), SD(i), Q] = get_avg_sd(V_ref(M))

% fprintf('-- IntVarP --');
% quantile(IntVarP(M), Q)

% tp - both the method and the rater found BGIDs
% tn - both the method and the rater did not find BGIDs
% ...
M_tp = M_vol & M_vol_ref;
M_fp = M_vol & ~M_vol_ref;
M_tn = ~M_vol & ~M_vol_ref;
M_fn = ~M_vol & M_vol_ref;

fprintf('-- FP volume --');
%quantile(V(M_fp), Q)
i=9; [Avg(i), SD(i), Q] = get_avg_sd(V(M_fp))
save_xls('FP', {'%s'}, Subjects(M_fp));

fprintf('-- FN volume --');
%quantile(V(M_fn), Q)
i=10; [Avg(i), SD(i), Q] = get_avg_sd(V_ref(M_fn))
save_xls('FN', {'%s'}, Subjects(M_fn));

fprintf('-- Confusion --');
Conf = [[sum(M_tp) sum(M_fp)];
        [sum(M_fn) sum(M_tn)]]
Sens = Conf(1, 1) / sum(Conf(:, 1));
Spec = Conf(2, 2) / sum(Conf(:, 2));

fprintf('-- Sens/Spec --');
[Sens Spec]

% Ba and Mba only with tp subjects
Subjects = Subjects(M);
figure;
Idx = blandAltmanPlot(V(M), V_ref(M))';
title(sprintf('N=%d', sum(M)));
% set(gca, 'XScale', 'log');
save_xls('BA_OLI', {'%s'}, Subjects(Idx));
set(gcf, 'color', 'white')

figure;
Idx = blandAltmanPlot(V(M), V_ref(M), J(M))';
title(sprintf('N=%d', sum(M)));
% set(gca, 'XScale', 'log');
save_xls('MBA_OLI', {'%s'}, Subjects(Idx));
set(gcf, 'color', 'white')

% save('forR', 'V', 'V_ref', 'J', '-v6');
% system('R CMD BATCH myqr.R');

function [Avg, SD, Q] = get_avg_sd(v)
Q = quantile(v(:), [.25 .5 .75]);
Avg = Q(2);
SD = (Q(3)-Q(1))/1.349;

