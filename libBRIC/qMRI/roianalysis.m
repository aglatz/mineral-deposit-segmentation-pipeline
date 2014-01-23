addpath('NIFTI/');
addpath('libBRIC/misc-matlab');

close all; clear all;

%Dir = '/home/aglatz/sf_V_DRIVE/Andreas_PhD/QMRI/Phantom/MnCl2/19676';
%S_roi = load_series(fullfile(Dir, 'Meta', 'DESPOT1_T1Map_roi'), []);
%S     = 1000./load_series(fullfile(Dir, 'R1_2', 'DESPOT1_T1Map'), []);
%S     = load_series(fullfile(Dir, '9', 'R2s'), []);
%Mat   = csvread(fullfile(Dir, 'Meta', 'Conc.csv'));

Dir = '/home/aglatz/sf_V_DRIVE/Andreas_PhD/QMRI/Phantom/MnCl2/19941';
S_roi = load_series(fullfile(Dir, 'Meta', 'DESPOT1HIFI_T1Map_roi'), []);
S = double(load_series(fullfile(Dir, '23', 'Phu_phuf0'), []));
Mat = csvread(fullfile(Dir, 'Meta', 'R1_map_conc.csv'));

% Dir = '/home/aglatz/sf_V_DRIVE/Andreas_PhD/QMRI/Phantom/MnCl2/19928';
% S_roi = load_series(fullfile(Dir, 'Meta', 'R1_map_roi'), []);
% S = double(load_series(fullfile(Dir, '3', 'R2'), []));
% Mat = csvread(fullfile(Dir, 'R1Dicom', 'Roi', 'R1_map_conc.csv'));

% figure;
% slice = round(size(S, 3)/2);
% S_tmp = S(:, :, slice);
% I_max = cast(quantile(double(S_tmp(:)), .8), class(S));
% plot_image_with_masks(S(:, :, slice), '', logical(S_roi(:, :, slice)), [], [], false, I_max);

Lab = Mat(:, 1);
N_Lab = length(Lab);
X = Mat(:, 2)';

idx = 1;
Slices = zeros(size(S, 3), 1);
for slice = 1:size(S, 3)
    I_tmp = S_roi(:, :, slice);
    Lab_tmp = unique(I_tmp(:));
    N = length(Lab_tmp) - 1;
    if N_Lab <= N
        Slices(idx) = slice;
        idx = idx + 1;
    end
end
Slices(idx:end) = [];
% Slices = 23:29; %7:16;
N_slices = length(Slices);

H1 = figure; hold on;
col = hsv(N_slices);
Res = zeros(N_slices, 2);
Y = zeros(N_slices, N_Lab);
Y_std = zeros(N_slices, N_Lab);
for slice_idx = 1:N_slices
	I = S(:, :, Slices(slice_idx));
    I_roi  = S_roi(:, :, Slices(slice_idx));
    for idx = 1:N_Lab
        M = I_roi == Lab(idx);
        I_tmp = I(M);
        I_tmp(isnan(I_tmp)) = [];
        if ~isempty(I_tmp)
            Y(slice_idx, idx) = median(I_tmp); %mloclogist(I_tmp);
            %Y(slice_idx, idx) = 1000./mloclogist(I_tmp);
            Y_std(slice_idx, idx) = iqr(I_tmp)./2; %mscalelogist(I_tmp);
        end
    end
    figure(H1);
    plot(X, Y(slice_idx, :), 'Color', col(slice_idx, :));
    %text(X, Y(slice_idx, :), num2str(ones(length(X), 1) .* slice_idx));
    errorbar(X, Y(slice_idx, :), Y_std(slice_idx, :), 'Color', col(slice_idx, :));
    
    Res(slice_idx, :) = polyfit(X, Y(slice_idx, :), 1);
end
xlabel('MnCl_2 concentration c in mMol');
ylabel('Mean ROI relaxivitiy R in s^{-1}');
Tmp_mean=median(Res);
Tmp_std=iqr(Res)./2;
text(0.06, 0, sprintf('R=((%0.1f\\pm%0.1f)mMol^{-1}c+(%0.1f\\pm%0.1f))s^{-1}', Tmp_mean(1), Tmp_std(1), Tmp_mean(2), Tmp_std(2)));
plot(X, polyval([Tmp_mean(1) Tmp_mean(2)], X), 'k', 'LineWidth', 3);


% set(gcf, 'color', 'w');
% export_fig([Name '.pdf'], '-a1',  '-q101');
