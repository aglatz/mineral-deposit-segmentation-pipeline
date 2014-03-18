addpath('NIFTI/');
addpath('libBRIC/misc-matlab');

close all; clear all;

Dir = '/home/aglatz/sf_V_DRIVE/Andreas_PhD/QMRI/Phantom/MnCl2/22218';
Rname = fullfile(Dir, 'irseepi', 'R1');
load([Rname '.mat']);
Y_ref = Y;
Y_ref_std = Y_std;
Fname = fullfile(Dir, 'R1_2fa', 'DESPOT1_T1Map');
load([Fname '.mat']);

Type='2';
figure;
%subplot(1,2,1);
hold on;
scatter(Y_ref, Y, 'filled');
xlabel(['\bf Reference relaxation rate R_{' Type ',ref} in s^{-1}']);
ylabel(['\bf Relaxation rate R_{' Type '} in s^{-1}']);
P = robustfit(Y_ref, Y);
Yest = polyval([P(2) P(1)], Y_ref);
plot(Y_ref, Yest, 'k'); %, 'LineWidth', 3);
MinV = min([Y_ref; Y]);
MaxV = max([Y_ref; Y]);
plot([MinV MaxV], [MinV MaxV], '--k');
if P(1) < 0
    Tmp = '-';
else
    Tmp = '+';
end 
text(MinV, MaxV, sprintf('R_{%s}=%0.2f R_{%s,ref} %s %0.2f', Type, P(2), Type, Tmp, abs(P(1))));
set(gcf, 'color', 'w');

figure
%subplot(1,2,2);
hold on;
Diff = (Y - Y_ref); %./Y_ref*100;
Avg =  (Y_ref + Y) / 2;
scatter(Avg, Diff, 'filled');
MeanDiff = median(Diff); mloclogist(Diff); median(Diff);
%CI_Diff = quantile(Diff, [.05 .95]);
SDDiff = iqr(Diff)/1.349; mscalelogist(Diff);
CI_Diff = [MeanDiff-SDDiff*1.95 MeanDiff+SDDiff*1.95];
MinAvg = min(Avg);
MaxAvg = max(Avg);
[P, stats] = robustfit(Avg, Diff);
ss = mean(abs(stats.resid))*sqrt(pi/2);
H(1) = plot([MinAvg MaxAvg], polyval([P(2) P(1)], [MinAvg MaxAvg]), 'k'); %[MeanDiff MeanDiff], 'k');
H(2) = plot([MinAvg MaxAvg], polyval([P(2) P(1)], [MinAvg MaxAvg])+stats.s*1.95, ':k'); %[CI_Diff(2) CI_Diff(2)], ':k');
legend(H, 'Bias', sprintf('95%% range (%0.2fs^{-1})', ss*1.95*2));
plot([MinAvg MaxAvg], polyval([P(2) P(1)], [MinAvg MaxAvg])-ss*1.95, ':k'); %[CI_Diff(1) CI_Diff(1)], ':k');
xlabel(['\bf (R_{' Type '}-R_{' Type ',ref})/2 in s^{-1}']);
ylabel(['\bf R_{' Type '}-R_{' Type ',ref} in s^{-1}']);
% off=0.02;
% text(MaxAvg*.8, MeanDiff+off, sprintf('Mean=%0.2fs^{-1}', MeanDiff));
% text(MaxAvg*.8, CI_Diff(1)+off, sprintf('Mean-1.95*SD=%0.2fs^{-1}', CI_Diff(1)));
% text(MaxAvg*.8, CI_Diff(2)+off, sprintf('Mean+1.95*SD=%0.2fs^{-1}', CI_Diff(2)));
set(gcf, 'color', 'w');
% Diff = [9, 6, -19, -8, -3, 7, 18, 42, 56, 58, 83, 69, 104, 103, 49, 21, -2, 91, -11];
% MeanDiff = mean(Diff); % median(Diff);
% SeDiff = std(Diff)/sqrt(length(Diff)); % iqr(Diff)/1.349/sqrt(length(Diff));
% dCI = tinv(.95, length(Diff)-1)*SeDiff;
% SeDiff1 = sqrt(3)*SeDiff;
% dCI1 = tinv(.95, length(Diff)-1)*SeDiff1;
% MinAvg = min(Avg);
% MaxAvg = max(Avg);
% plot([MinAvg MaxAvg], [MeanDiff MeanDiff], 'k');
% plot([MinAvg MaxAvg], [MeanDiff+dCI MeanDiff+dCI], '--k');
% plot([MinAvg MaxAvg], [MeanDiff-dCI MeanDiff-dCI], '--k');
% plot([MinAvg MaxAvg], [MeanDiff+dCI1 MeanDiff+dCI1], 'k');
% plot([MinAvg MaxAvg], [MeanDiff-dCI1 MeanDiff-dCI1], 'k');

% 
% csvwrite([Fname '.csv'], [X Y Y_std]);
% save([Fname '.mat'], 'C', 'X', 'Y', 'Y_std');
% Mat = cell2mat(C); %(2:end-1, :));
% [p, tbl, stats]=kruskalwallis(Mat(:, 1), Mat(:, 2), 'off');
% figure; Res = multcompare(stats);
% Tmp = sign(Res(:, 3)) + sign(Res(:, 5));
% tbl
% Res(Tmp < 2 & Tmp > -2, 1:2)-1


% set(gcf, 'color', 'w');
% export_fig([Name '.pdf'], '-a1',  '-q101');
