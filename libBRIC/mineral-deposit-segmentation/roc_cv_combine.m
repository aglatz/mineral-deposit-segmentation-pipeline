addpath ../../libBRIC/misc-matlab/

close all; %clear all; clc

SubjectFile = 'subjects_98';
Pfx = 'ad_t1w_mineral_stdint'; % varies the T1w hypointensity factor
AdaptiveFlag = true;
RoiLabelTable = {[13, 11, 12, 14]};
N_cpus = 6;
Thr = 0:0.1:1.5;
N_thr = length(Thr);
N_cv = 10;
out = struct;

H2 = figure;
Col = hsv(N_cv);
J_mean = NaN(N_cv, N_thr);
for idx_dir = 1:N_cv
    if idx_dir == 9
        fprintf('');
    end
    
    for idx_thr = 1:N_thr
        [~, J] = load_matdata(fullfile(num2str(idx_dir), ['hu_' num2str(Thr(idx_thr)) '_' Pfx '.mat']));
        J_mean(idx_dir, idx_thr) = quantile(J, .5);
    end    

	idx = find(J_mean(idx_dir, :) == max(J_mean(idx_dir, :)));
    [idx, len] = find_largemax(idx, '');
%     if mod(len, 2) == 0
%         idx = [idx-1, idx];
%     end

    figure(H2); hold on; plot(Thr, J_mean(idx_dir, :), 'color', Col(idx_dir, :));
    out(idx_dir).thr = mean(Thr(idx)); 
    out(idx_dir).J = (J_mean(idx_dir, idx)-J_mean(idx_dir, 1))./J_mean(idx_dir, 1);
    scatter(out(idx_dir).thr, (J_mean(idx_dir, idx)), 'k', 'filled');
    
%     SubjectFile_traintest = [SubjectFile '_' Pfx '_' num2str(idx_dir)];
%     [out(idx_dir).Ret, out(idx_dir).Subjects] = ...
%                     segment_us_mp([SubjectFile_traintest '.xls'], ...
%                                   RoiLabelTable, N_cpus, ...
%                                   'ThreshFactor', [1 0], ...
%                                   'AdaptiveFlag', AdaptiveFlag, ...
%                                   'SaveMaskFlag', true, ...
%                                   'ReportName', 'class', ...
%                                   'N_gre', 1, ...
%                                   'CNR_thr', 0, ...
%                                   'phypo_thr', 0.1, ...
%                                   'intvar_thr', out(idx_dir).thr);

end

% MatFileName = [SubjectFile '_' Pfx '.mat'];
% save(MatFileName, 'out');

figure(H2);
axis([0 1.6 0.295 0.705]);
Avg=median([out.thr]);
SD=iqr([out.thr])./ 1.349;
plot([Avg-SD Avg-SD], [0.295 0.705], '--k', 'LineWidth', 1);
plot([Avg+SD Avg+SD], [0.295 0.705], '--k', 'LineWidth', 1);
plot([Avg Avg], [0.295 0.705], 'k', 'LineWidth', 2);
%plot(Thr, quantile(J_mean, 1-.1587), '--k', 'LineWidth', 2);
%plot(Thr, quantile(J_mean, .1587), '--k', 'LineWidth', 2);
xlabel('\bf Local variance filter parameter');
ylabel('\bf Jaccard index');
set(gcf, 'color', 'white');

