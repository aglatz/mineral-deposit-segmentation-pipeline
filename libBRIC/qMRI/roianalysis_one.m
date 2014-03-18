addpath('NIFTI/');
addpath('LIBRA/');
addpath('libBRIC/misc-matlab');

close all; clear all; clc

Dir = '/home/aglatz/sf_V_DRIVE/Andreas_PhD/QMRI/Phantom/MnCl2/22218';
S_roi = load_series(fullfile(Dir, 'Meta', 'R2s_roi'), []);
Fname = fullfile(Dir, '16', 'R2s');
S = double(load_series(Fname, []));
Mat = csvread(fullfile(Dir, 'Meta', 'R1_map_conc.csv'));

% Dir = '/home/aglatz/sf_V_DRIVE/Andreas_PhD/QMRI/Subjects/19860';
% S_roi = load_series(fullfile(Dir, 'Roi', 'T1w_R1HIFI_roi'), []);
% Fname = fullfile(Dir, 'R1HIFI', 'DESPOT1HIFI_T1Map');
% S = 1000./double(load_series(Fname, []));
% Tmp = 0:7;
% Mat = [2.^Tmp', Tmp'];

% figure;
% slice = round(size(S, 3)/2);
% S_tmp = S(:, :, slice);
% I_max = cast(quantile(double(S_tmp(:)), .8), class(S));
% plot_image_with_masks(S(:, :, slice), '', logical(S_roi(:, :, slice)), [], [], false, I_max);

cal = 3;
Type = '2*';
off = 0.1;
switch cal
    case 1,
        S = 1.30*S-0.09;
    case 2,
        S = 1.02*S+0.24;
    case 3,
        S = 0.97*S+0.86;
    otherwise,
end

Lab = Mat(:, 1);
N_Lab = length(Lab);
X = Mat(:, 2);

Y = zeros(N_Lab, 1);
Y_std = zeros(N_Lab, 1);
C = cell(N_Lab, 2);
for idx = 1:N_Lab
    SM = S_roi == Lab(idx);
    I_tmp = S(SM);
    I_tmp(isnan(I_tmp)) = [];
    if ~isempty(I_tmp)
        Y(idx) = median(I_tmp); %mloclogist(I_tmp); %median(I_tmp);
        Y_std(idx) = iqr(I_tmp)/1.349; %mscalelogist(I_tmp); %iqr(I_tmp)./2;
        C{idx, 1} = I_tmp;
        C{idx, 2} = ones(size(I_tmp))*idx;
    end
end

figure; hold on;
errorbar(X, Y, Y_std);
xlabel('\bf MnCl_2 concentration c in mmol/l');
ylabel(['\bf Relaxivitiy rate R_{' Type '} in s^{-1}']);
Xmin = min(X);
Xmax = max(X);
[P, stats] = robustfit(X, Y);
Yest = polyval([P(2) P(1)], X);
plot(X, Yest, 'k', 'LineWidth', 3);
if P(1) < 0
    Tmp = '-';
else
    Tmp = '+';
end 
text(Xmin, max(Y), sprintf(['R_{' Type '}=%0.2f s^{-1}/mmol/l c %s %0.2f s^{-1}'], P(2), Tmp, abs(P(1))));
set(gcf, 'color', 'w');

figure; hold on;
e = stats.rstud;
qqplot(e);

figure; hold on;
plot(X, e);
ME = mloclogist(e); %median(e);
plot([Xmin Xmax], [ME ME], 'k');
SD = mscalelogist(e); %iqr(e)/1.349;
plot([Xmin Xmax], [ME+SD ME+SD], '--k');
plot([Xmin Xmax], [ME-SD ME-SD], '--k');
%t025=tinv(.025, length(e)-1);
t975=tinv(.975, length(e)-1);
plot([Xmin Xmax], [ME+t975 ME+t975], ':k');
plot([Xmin Xmax], [ME-t975 ME-t975], ':k');
xlabel('\bf MnCl_2 concentration c in mmol/l');
ylabel('\bf Standardised residuals');
% r2 = corr(Y, P(1) + P(2)*X)^2;
% sse = stats.dfe * stats.robust_s^2;
% yhat = P(1) + P(2)*X;
% ssr = norm(yhat-mean(yhat))^2;
% r2 = 1 - sse / (sse + ssr);
Yhat = P(1) + P(2)*X;
nrms = 1-norm(Y-Yhat)/norm(Y-mean(Yhat));
title(sprintf('%0.2f;%0.2f;nrms=%0.2f', ME, SD, nrms));
text(Xmax*.8, ME+off, sprintf('Mean=%0.2f s^{-1}', ME));
text(Xmax*.8, ME+SD+off, sprintf('Mean+SD=%0.2f s^{-1}', ME+SD));
text(Xmax*.8, ME-SD+off, sprintf('Mean-SD=%0.2f s^{-1}', ME-SD));
set(gcf, 'color', 'w');

csvwrite([Fname '.csv'], [X Y Y_std]);
save([Fname '.mat'], 'C', 'X', 'Y', 'Y_std');
Mat = cell2mat(C); %(2:end-1, :));
[p, tbl, stats]=kruskalwallis(Mat(:, 1), Mat(:, 2), 'off');
figure; Res = multcompare(stats);
Tmp = sign(Res(:, 3)) + sign(Res(:, 5));
tbl
Res(Tmp < 2 & Tmp > -2, 1:2)-1



% set(gcf, 'color', 'w');
% export_fig([Name '.pdf'], '-a1',  '-q101');
