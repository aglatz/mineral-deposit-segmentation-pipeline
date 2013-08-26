close all; clear all; clc

Name = 'subjects_98';
Pfx = 'ad_t1w0.975';
N_cpus = 6;
AdaptiveFlag = true;
RoiLabelTable = {[13 11 12 14]};
% [Ret, Subjects, Over] = segment_us_mp([Name '.xls'], ...
%                                       RoiLabelTable, N_cpus, ...
%                                       'ThreshFactor', [1 0], ...
%                                       'AdaptiveFlag', AdaptiveFlag, ...
%                                       'SaveMaskFlag', true, ...
%                                       'ReportName', 'class', ...
%                                       'IntvarP', 0.4);
% save([Name '_' Pfx '.mat']);

[Subjects, J, D, V, V_ref, IntVarP] = load_matdata([Name '_' Pfx '.mat']);

fprintf('-- Jaccard (All) --');
Q = [.25 .5 .75];
quantile(J, Q)

fprintf('-- Dice (All) --');
quantile(D, Q)

fprintf('-- Volume (All) --');
quantile(V, Q)

fprintf('-- Reference volume (All) --');
quantile(V_ref, Q)

fprintf('-- IntVarP (All) --');
quantile(IntVarP, Q)

M_vol = V > 0;
M_vol_ref = V_ref > 0;
M = M_vol & M_vol_ref;

fprintf('-- Jaccard --');
Q = [.25 .5 .75];
quantile(J(M), Q)

fprintf('-- Dice --');
quantile(D(M), Q)

fprintf('-- Volume --');
quantile(V(M), Q)

fprintf('-- Reference volume --');
quantile(V_ref(M), Q)

fprintf('-- IntVarP --');
quantile(IntVarP(M), Q)

% tp - both the method and the rater found BGIDs
% tn - both the method and the rater did not find BGIDs
% ...
M_tp = M_vol & M_vol_ref;
M_fp = M_vol & ~M_vol_ref;
M_tn = ~M_vol & ~M_vol_ref;
M_fn = ~M_vol & M_vol_ref;

fprintf('-- FP volume --');
quantile(V(M_fp), Q)
save_xls('FP', Subjects(M_fp));

fprintf('-- FN volume --');
quantile(V(M_fn), Q)
save_xls('FN', Subjects(M_fn));

fprintf('-- Confusion --');
Conf = [[sum(M_tp) sum(M_fp)];
        [sum(M_fn) sum(M_tn)]]
Sens = Conf(1, 1) / sum(Conf(:, 1));
Spec = Conf(2, 2) / sum(Conf(:, 2));

fprintf('-- Sens/Spec --');
[Sens Spec]

% Ba and Mba only with tp subjects
Subjects = Subjects(M_tp);
figure;
Idx = blandAltmanPlot(V(M_tp), V_ref(M_tp))';
set(gca, 'XScale', 'log');
save_xls('BA_OLI', Subjects(Idx));
set(gcf, 'color', 'white')

figure;
Idx = blandAltmanPlot(V(M_tp), V_ref(M_tp), J(M_tp))';
set(gca, 'XScale', 'log');
save_xls('MBA_OLI', Subjects(Idx));
set(gcf, 'color', 'white')

% save('forR', 'V', 'V_ref', 'J', '-v6');
% system('R CMD BATCH myqr.R');

