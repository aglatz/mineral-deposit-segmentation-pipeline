addpath('NIFTI/');
addpath('libBRIC/misc-matlab');

close all; clear all;

%Dir = '/home/aglatz/sf_V_DRIVE/Andreas_PhD/QMRI/Phantom/MnCl2/19676';
%S_roi = load_series(fullfile(Dir, 'Meta', 'DESPOT1_T1Map_roi'), []);
%S     = 1000./load_series(fullfile(Dir, 'R1_2', 'DESPOT1_T1Map'), []);
%S     = load_series(fullfile(Dir, '9', 'R2s'), []);
%Mat   = csvread(fullfile(Dir, 'Meta', 'Conc.csv'));

Dir = '/home/aglatz/sf_V_DRIVE/Andreas_PhD/QMRI/Phantom/MnCl2/19906';
S_roi = load_series(fullfile(Dir, 'Meta', 'DESPOT1_T1Map_roi'), []);
Fname = fullfile(Dir, '20', 'R2s');
S = double(load_series(Fname, []));
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
X = Mat(:, 2);

Y = zeros(N_Lab, 1);
Y_std = zeros(N_Lab, 1);
for idx = 1:N_Lab
    SM = S_roi == Lab(idx);
    I_tmp = S(SM);
    I_tmp(isnan(I_tmp)) = [];
    if ~isempty(I_tmp)
        Y(idx) = median(I_tmp); %mloclogist(I_tmp); %median(I_tmp);
        Y_std(idx) = iqr(I_tmp)/1.349; %mscalelogist(I_tmp); %iqr(I_tmp)./2;
    end
end

figure; hold on;
errorbar(X, Y, Y_std);
xlabel('\bf MnCl_2 concentration in mMol');
ylabel('\bf Mean ROI relaxivitiy rate in s^{-1}');

P = robustfit(X, Y);
Yest = polyval([P(2) P(1)], X);
text(0.04, max(Yest), sprintf('P=(%0.2fmMol^{-1}s^{-1};%0.2fs^{-1})', P(2), P(1)));
plot(X, Yest, 'k', 'LineWidth', 3);

csvwrite([Fname '.csv'], [X Y Y_std]);


% set(gcf, 'color', 'w');
% export_fig([Name '.pdf'], '-a1',  '-q101');
