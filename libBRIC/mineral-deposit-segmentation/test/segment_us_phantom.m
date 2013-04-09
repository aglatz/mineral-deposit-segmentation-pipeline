% First include our libraries
addpath('../../mineral-deposit-segmentation');
addpath('../../misc-matlab/');

% Include external libraries
addpath('/home/aglatz/tmp/mineral/NIFTI/');
addpath('/home/aglatz/tmp/mineral/LIBRA/');

close all; clear all; clc

% File for summarizing classification statistics
OutFile = 'class_us';
Subject = 'data';
Out_name = [Subject '/' OutFile];
save_ps_figure(Out_name, []); % Delete previous file

% Read ROI mask
Name_roi_red = 'RO_mask';
Name_roi = [Subject '/' Name_roi_red];
S_roi = load_series(Name_roi, []);
Roi = roi_init(S_roi);
S_roi = load_series(Name_roi, roi_nifti_sliceno(Roi, []));

% Read reference mask if exists
try
    Fe_name = 'FE_roi_mask';
    S_ref = load_series([Subject '/' Fe_name], roi_nifti_sliceno(Roi, []));
catch
    S_ref = [];
end

% Read co-registered T1- and T2*-weighted volumes
S_gre = double(load_series([Subject '/GRE_restore'], roi_nifti_sliceno(Roi, [])));
S_t1w = double(load_series([Subject '/T1W_restore'], roi_nifti_sliceno(Roi, [])));

% Eliminate zeros
% SM = logical(S_roi);
% M = S_t1w(SM) <= eps | S_gre(SM) <= eps;
% SM_tmp = SM;
% SM_tmp(SM) = M;
% S_roi(SM_tmp) = zeros(sum(SM_tmp(:)), 1);
    
% Get labels
Lab = unique(S_roi)';
Lab = Lab(2:end); % Exclude background
N_lab = length(Lab);

% Init output variables
S_out_all = zeros(size(S_roi), class(S_roi));
S_nontis_all = zeros(size(S_roi), class(S_roi));
S_ntis_all = zeros(size(S_roi), class(S_roi));
I_means_est = zeros(N_lab, 2);
I_ntis_means_all = zeros(N_lab, 2);

% MAIN loop: go trough all the labels
for idx_lab = 1:N_lab
    % Estimate the normal-appearing tissue means (agarose)
    SM_tmp = S_roi == Lab(idx_lab);
	Mat = [S_gre(SM_tmp) S_t1w(SM_tmp)];
	[~, I_means_est(idx_lab, :)] = pcomp_find(Mat);
        
    % Segmentation
    [S_out, S_nontis, S_ntis, I_ntis_means_all(idx_lab, :), I_gre_thr, I_gre_thr_ref] = ...
        segment_us(S_gre, S_t1w, S_roi, S_ref, Lab(idx_lab), I_means_est(idx_lab, :), ...
                   [], 'Out_name', Out_name); %, 'P_thr', [0.9579 0.4314]);
    S_out_all = S_out_all + S_out;
    S_nontis_all = S_nontis_all + S_nontis;
    S_ntis_all = S_ntis_all + S_ntis;
end
    
% Summary plot
H = figure; %create_ps_figure; %figure;
Col = [[0, 0, 0]; hsv(3)];
scatter(S_gre(logical(S_roi)), S_t1w(logical(S_roi)), 10, Col(1, :));
hold on;
scatter(I_ntis_means_all(:, 1), I_ntis_means_all(:, 2), 20, 'r');
scatter(I_means_est(:, 1), I_means_est(:, 2), 20, 'g');
SM_out = logical(S_out_all);
Mat = [S_gre(SM_out) S_t1w(SM_out)];
scatter(Mat(:, 1), Mat(:, 2), 10, Col(2, :));
if ~isempty(S_ref)
    SM_ref = logical(S_ref);
    scatter(S_gre(SM_ref & SM_out), S_t1w(SM_ref & SM_out), 10, Col(3, :));
    scatter(S_gre(SM_ref & ~SM_out), S_t1w(SM_ref & ~SM_out), 10, Col(4, :));
end
title(sprintf('Sum: %d', sum(SM_out(:))));
axis equal;
save_ps_figure(Out_name, H);

% Save results
Out_name = 'Oli_usnew_mask';
save_series(Name_roi, [Subject '/' Out_name], S_out_all, roi_nifti_sliceno(Roi, []));
save_series(Name_roi, [Subject '/NonTis_mask'], S_nontis_all, roi_nifti_sliceno(Roi, []));
save_series(Name_roi, [Subject '/NormTis_mask'], S_ntis_all, roi_nifti_sliceno(Roi, []));

if ~isempty(S_ref)
    Ret = validate(Subject, Out_name, Fe_name, [], Name_roi_red, 1);
    save([Subject '/Ret.mat'], 'Ret');
end

