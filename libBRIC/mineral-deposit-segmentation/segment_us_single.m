function [Ret] = segment_us_single(Subject, RoiLabelTable, ReportName, ...
                                   InterpFactor, ThreshFactor, ...
                                   AdaptiveFlag, varargin)
% This function expects three files in each subject directory:
% - RO_mask: the masks with the regions of interests
% - GRE_brain_restore: the GRE volume which was brain extracted
%                      (not strictly necessary) and bias-field corrected.
% - T1W_brain_restore: the T1W volume which was brain extracted
%                      (not strictly necessary) and bias-field corrected
% It can generate following masks for the given subject unless the optional
% argument 'SaveMaskFlag' is 'false':
% - NonTis_mask: A mask selecting what appears to be not normal tissue
% - NormTis_mask: A mask selecting the normal-appearing tissue
% - T2swHypo_mask: A mask selecting T2sw hypointensities (subset of
%                  NonTis_mask)
% It can save the intermediate plots to a file if 'ReportName' is
% a valide file name.
% 
% INPUTS: Subject - the subject directory containing at least RO_mask,
%                   GRE_brain_restore, T1W_brain_restore
%         RoiLabelTable - A cell vector whose elements are label vectors.
%                         The first entry of the label vector are the
%                         labels of the ROIs that should be used to
%                         estimate the thresholds. All other ROI labels of
%                         the same vector should be segmented with the same
%                         threshold.
%         InterpFactor - Controls FFT interpolation for greyscale volumes
%                        and nearest neighbour interpolation for mask
%                        volumes. A 'InterpFactor'=1 means no
%                        interpolation.
%         ThreshFactor - T2*w threshold refinement parameters (see
%                        segment_us())
%         AdaptiveFlag - If set to 'true' the T2*w is refined with
%                        an adaptive method
% OPTIONAL INPUTS: SaveMaskFlag - Controls the saving of the generated
%                                 masks. Default is 'true'.
% OUTPUT: The structure from validate_raw() plus additional element
%         which is the T2*w threshold.
%

% Process optional args
N_vain = length(varargin);
SaveMaskFlag = true;
for idx_vain = 1:N_vain
    arg_in = varargin{idx_vain};
    if ischar(arg_in) || isscalar(arg_in)
        switch arg_in
            case {'SaveMaskFlag'}
                SaveMaskFlag = varargin{idx_vain+1};
        end
    end
end

% Delete previous version of report file
if ~isempty(ReportName)
    ReportFile = [Subject '/' ReportName];
    save_ps_figure(ReportFile, []); % Deletes previous file
else
    ReportFile = [];
end

% Read ROI
RoiName = 'RO_mask';
RoiFile = [Subject '/' RoiName];
S_roi = load_series(RoiFile, []);
Roi = roi_init(S_roi);
S_roi = load_series_interp(RoiFile, roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);

% Pool signal intensities from corresponding left and right hemisphere structures
S_roi(S_roi == 50) = 11;
S_roi(S_roi == 51) = 12;
S_roi(S_roi == 52) = 13;
S_roi(S_roi == 55) = 14;

% Read Reference
try
    FeName = 'FE_roi_mask';
    S_ref = load_series_interp([Subject '/' FeName], roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);
catch
    S_ref = [];
end

% Read T2sw/T1w volumes
T2swName = 'GRE_brain_restore';
S_gre = double(load_series_interp([Subject '/' T2swName], roi_nifti_sliceno(Roi, []), 'fft', InterpFactor));
T1wName = 'T1W_brain_restore';
S_t1w = double(load_series_interp([Subject '/' T1wName], roi_nifti_sliceno(Roi, []), 'fft', InterpFactor));

% Initialize output variables and start segmentation
N_iter = length(RoiLabelTable);
S_out_all = zeros(size(S_roi), class(S_roi));
S_out_hypo_all = zeros(size(S_roi), class(S_roi));
S_out_hyper_all = zeros(size(S_roi), class(S_roi));
S_nontis_all = zeros(size(S_roi), class(S_roi));
S_ntis_all = zeros(size(S_roi), class(S_roi));
I_ntis_means = cell(N_iter, 1);
I_gre_thr = zeros(N_iter, 1);
for idx_iter = 1:N_iter
    Labs = RoiLabelTable{idx_iter};
%     N_labs = length(Labs);
%
%         % Heuristic to estimate the normal-appearing tissue means - works
%         % for GRE and T1W of LBC1936 protocol
%         I_ntis_means_all = NaN(N_labs, 2);
%         for idx_lab = 1:N_labs
%             SM_tmp = S_roi == Labs(idx_lab);
%             Mat = [S_gre(SM_tmp) S_t1w(SM_tmp)];
%             [~, I_ntis_means_all(idx_lab, :)] = pcomp_find(Mat);
%         end
%         Tmp = mcdregres(I_ntis_means_all(:, 2), I_ntis_means_all(:, 1), 'plots', 0);
%         Est_param = [Tmp.slope Tmp.int];
%         I_ntis_means_all_est = I_ntis_means_all;
%         I_ntis_means_all_est(:, 1) = polyval(Est_param, I_ntis_means_all(:, 2));

    % Segmentation
    Ret = segment_us(S_gre, S_t1w, S_roi, S_ref, Labs, [], [], AdaptiveFlag, ...
                     'Out_name', ReportFile, 'P_thr', ThreshFactor);

    % Store results
    S_out_all = S_out_all + Ret.S_out;
    S_out_hypo_all = S_out_hypo_all + Ret.S_out_hypo;
    S_out_hyper_all = S_out_hyper_all + Ret.S_out_hyper;
    S_nontis_all = S_nontis_all + Ret.S_nontis;
    S_ntis_all = S_ntis_all + Ret.S_ntis;
    I_ntis_means{idx_iter} = Ret.I_ntis_means;
    % Save for later...
    I_gre_thr(idx_iter) = Ret.I_gre_thr;
end

% Summary plot
H = figure; %create_ps_figure;
Col = [[0, 0, 0]; hsv(3)];
scatter(S_gre(logical(S_roi)), S_t1w(logical(S_roi)), 10, Col(1, :));
hold on;
for idx_iter = 1:N_iter
    I_ntis_means_iter = I_ntis_means{idx_iter};
    scatter(I_ntis_means_iter(:, 1), I_ntis_means_iter(:, 2), 20, 'r');
end
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
save_ps_figure(ReportFile, H);

% Change resolution of masks
T2swHypoMaskName = 'T2swHypo_mask';
S_out_all = interp_series([Subject '/' RoiName], S_out_all, ...
                roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);
S_out_hypo_all = interp_series([Subject '/' RoiName], S_out_hypo_all, ...
                roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);
S_out_hyper_all = interp_series([Subject '/' RoiName], S_out_hyper_all, ...
                roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);
S_ntis_all = interp_series([Subject '/' RoiName], S_ntis_all, ...
                roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);                
[SM_oli] = morph_filter(logical(S_out_all), S_gre, S_roi, logical(S_ntis_all), logical(S_out_hypo_all));
S_out_all(~SM_oli) = 0;
S_out_hypo_all(~SM_oli) = 0;
S_out_hyper_all(~SM_oli) = 0;

if SaveMaskFlag
    % Save masks
    save_series([Subject '/' RoiName], [Subject '/' T2swHypoMaskName], ...
                S_out_all, roi_nifti_sliceno(Roi, []));
    save_series([Subject '/' RoiName], [Subject '/T2swHypoT1wHypo_mask'], ...
                S_out_hypo_all, roi_nifti_sliceno(Roi, []));
    save_series([Subject '/' RoiName], [Subject '/T2swHypoT1wHyper_mask'], ...
                S_out_hyper_all, roi_nifti_sliceno(Roi, []));
    % Save non tissue mask
    save_series_interp([Subject '/' RoiName], [Subject '/NonTis_mask'], ...
                       S_nontis_all, roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);
    % Save norm tissue mask
    save_series([Subject '/' RoiName], [Subject '/NormTis_mask'], ...
                S_ntis_all, roi_nifti_sliceno(Roi, []));
end

% Validate
NII = load_series([Subject '/' RoiName], 0);
F = NII.hdr.dime.pixdim(2:4);
if ~isempty(S_ref)
    Ret = validate_raw(S_out_all, S_ref, [], F);
else
    % Fake validation just to get the volume
    Ret = validate_raw(S_out_all, S_out_all, [], F);
end
Ret.I_gre_thr = I_gre_thr; % Also add T2*w threshold