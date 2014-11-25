addpath('NIFTI/');
addpath('libBRIC/misc-matlab');

close all; clear all;

Dir = '/home/aglatz/sf_V_DRIVE/Andreas_PhD/QMRI/Phantom/MnCl2/22218';
Rname = fullfile(Dir, 'irseepi', 'R1'); %'se', 'R2');
load([Rname '.mat']);
Y_ref = Y;
Y_ref_std = Y_std;
% Fname = fullfile(Dir, 'R1_2fa', 'DESPOT1_T1Map'); %'17', 'R2'); %'16', 'R2s'); %
% load([Fname '.mat']);
Dir = '/home/aglatz/sf_V_DRIVE/Andreas_PhD/QMRI/Phantom/MnCl2/23025/';
Rname = fullfile(Dir, '8_9_10_11_12_13', 'DESPOT1HIFI_T1Map'); %'se', 'R2');
load([Rname '.mat']);

Type='1';
figure;
%subplot(1,2,1);
hold on;
scatter(Y_ref, Y, 'filled');
xlabel(['\bf R' Type '_{ref} in 1/s']);
ylabel(['\bf R' Type ' in 1/s']);
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
text(MinV, MaxV, sprintf(['\\bf R' Type '=%0.2f R' Type '_{ref} %s %0.2f 1/s'], P(2), Tmp, abs(P(1))));
set(gcf, 'color', 'w');
set(gca, 'FontSize', 12);
set(findall(gcf, 'type', 'text'), 'FontSize', 12);
export_fig([Fname '_equality.eps']);

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
H(3) = plot([MinAvg MaxAvg], [0 0], '--k');
legend(H(1:2), 'Bias', sprintf('95%% range (%0.2f 1/s)', ss*1.95*2));
plot([MinAvg MaxAvg], polyval([P(2) P(1)], [MinAvg MaxAvg])-ss*1.95, ':k'); %[CI_Diff(1) CI_Diff(1)], ':k');
xlabel(['\bf (R' Type '+R' Type '_{ref})/2 in 1/s']);
ylabel(['\bf R' Type '-R' Type '_{ref} in 1/s']);
set(gcf, 'color', 'w');
set(gca, 'FontSize', 12);
set(findall(gcf, 'type', 'text'), 'FontSize', 12);
export_fig([Fname '_ba.eps']);
