addpath ../../libBRIC/misc-matlab/

close all; clear all; clc

SubjectFile = 'subjects_98';
Pfx = 'fi_t1w0.975';
AdaptiveFlag = false;
RoiLabelTable = {[13, 11, 12, 14]};
N_cpus = 6;
Thr = (0:0.1:0.9)';
N_cv = 10;
out = struct;

H1 = figure;
H2 = figure;
Col = hsv(N_cv);
for idx_dir = 1:N_cv
    N_thr = length(Thr);
    J_mean = zeros(N_thr, 1);
    Sens = zeros(N_thr, 1);
    Spec = zeros(N_thr, 1);
    
    for idx_thr = 1:N_thr
        [J, ~, V, V_ref, ~, ~, ~, ~, Sens(idx_thr), Spec(idx_thr)] = ...
            load_matdata([num2str(idx_dir) '/hu_' num2str(Thr(idx_thr)) '_' Pfx '.mat'], 0);
        J_mean(idx_thr) = quantile(J, .5);
    end    

	idx = find(J_mean == max(J_mean));
    [idx, len] = find_largemax(idx, '');
    if mod(len, 2) == 0
        idx = [idx-1, idx];
    end
    figure(H1); hold on; plot(1-Spec, Sens, 'color', Col(idx_dir, :));
    x = 1 - mean(Spec(idx));
    y = mean(Sens(idx));
    scatter(x, y, 'k', 'filled');

    figure(H2); hold on; plot(Thr, J_mean, 'color', Col(idx_dir, :));
    out(idx_dir).thr = mean(Thr(idx));
    scatter(out(idx_dir).thr, mean(J_mean(idx)), 'k', 'filled');
    
%     SubjectFile_traintest = [SubjectFile '_' Pfx '_' num2str(idx_dir)];
%     [out_cc(idx_dir).Ret, out_cc(idx_dir).Subjects] = ...
%                     segment_us_mp([SubjectFile_traintest '.xls'], ...
%                                   RoiLabelTable, N_cpus, ...
%                                   'ThreshFactor', [1 0], ...
%                                   'ReportName', 'class', ...
%                                   'AdaptiveFlag', AdaptiveFlag, ...
%                                   'SaveMaskFlag', true, ...
%                                   'IntvarP', out_cc(idx_dir).thr);
end

figure(H1);
axis([0 1 0 1]);
plot([0 1], [0 1], '--k');
xlabel('\bf 1-Specificity');
ylabel('\bf Sensitivity');

figure(H2);
xlabel('\bf Cumulative probability');
ylabel('\bf Jaccard index');
