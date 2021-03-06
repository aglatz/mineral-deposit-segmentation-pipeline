function [Ret, CC] = segment_us_single(Subject, RoiLabelTable, ReportName, ...
                                       ThreshFactor, AdaptiveFlag, N_gre, ...
                                       CNR_thr, intvar_thr, varargin)
% This function expects three files in each subject directory:
% - RO_mask: the masks with the regions of interests
% - GRE_brain_restore: the GRE volume
% - T1W_brain_restore: the T1W volume
% - FE_roi_mask (optional): a reference mask for T2*w hypointensities
% It can generate following masks for the given subject unless the optional
% argument 'SaveMaskFlag' is 'false':
% - NonTis_mask: A mask selecting what appears to be not normal tissue
% - NormTis_mask: A mask selecting the normal-appearing tissue
% - T2swHypo_mask: A mask selecting T2sw hypointensities (subset of
%                  NonTis_mask)
% It can save the intermediate plots to a file if 'ReportName' is valid.
% 
% INPUTS: Subject - the subject directory containing at least RO_mask,
%                   GRE_brain_restore, T1W_brain_restore
%         RoiLabelTable - A cell vector whose elements are label vectors.
%                         The first entry of the label vector are the
%                         labels of the ROIs that should be used to
%                         estimate the thresholds. All other ROI labels of
%                         the same vector should be segmented with the same
%                         threshold.
%         ReportName - A valid file name or []; If file name is valid the
%                      generated plots are saved.
%         ThreshFactor - T2*w threshold refinement parameters (see
%                        segment_us())
%         AdaptiveFlag - If set to 'true' the T2*w is refined with
%                        an adaptive method
%         N_gre - Estimated noise of T2*w magnitude volume
%         CNR_thr - contrast-to-noise ratio threshold
%         intvar_thr - T2*w intensity variance threshold
% OPTIONAL INPUTS: SaveMaskFlag - Controls the saving of the generated
%                                 masks. Default is 'true'.
% OUTPUT: The structure from validate_raw() plus additional element
%         which is the T2*w threshold.
%
% EXAMPLE:
%   addpath('LIBRA/') % e.g. symbolic link to LIBRA toolbox
%   addpath('NIFTI/') % e.g. symbolic link to NIFTI toolbox
%   Subject = '/home/aglatz/tmp/mineral/2/13779';
%   RoiLabelTable = {[13 11 12 14]};
%   ThreshFactor = [1 0];
%   AdaptiveFlag = true;
%   ReportName = 'class';
%   IntvarP = 0.45;
%   InterpFactor = 1;
%   segment_us_single(Subject, RoiLabelTable, ReportName, ...
%                     InterpFactor, ThreshFactor, ...
%                     AdaptiveFlag, IntvarP);
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
% Subject = '~/tmp/mineral/4/14073';
if ~isempty(ReportName)
    ReportFile = fullfile(Subject, ReportName);
    save_ps_figure(ReportFile, []); % Deletes previous file
else
    ReportFile = [];
end

% Read ROI
RoiName = 'RO_mask';
RoiFile = fullfile(Subject, RoiName);
S_roi = load_series(RoiFile, []);
Roi = roi_init(S_roi);
S_roi = load_series(RoiFile, roi_nifti_sliceno(Roi, []));

% Pool signal intensities from corresponding left and right hemisphere structures
% S_roi(S_roi == 50) = 11;
% S_roi(S_roi == 51) = 12;
% S_roi(S_roi == 52) = 13;
% S_roi(S_roi == 55) = 14;

% Read Reference
try
    FeName = 'FE_roi_mask';
    S_ref_orig = load_series(fullfile(Subject, FeName), roi_nifti_sliceno(Roi, []));
%     S_ref_orig(S_ref_orig == 50) = 11;
%     S_ref_orig(S_ref_orig == 51) = 12;
%     S_ref_orig(S_ref_orig == 52) = 13;
%     S_ref_orig(S_ref_orig == 55) = 14;
    S_ref = zeros(size(S_ref_orig), class(S_ref_orig));
    N_iter = length(RoiLabelTable);
    for idx_iter = 1:N_iter
        Labs = RoiLabelTable{idx_iter};
        for lab = Labs
            S_ref(S_ref_orig == lab) = lab;
        end
    end
catch
    S_ref = [];
end

% Read T2sw/T1w volumes
T2swName = 'GRE_brain_restore';
S_gre = double(load_series(fullfile(Subject, T2swName), roi_nifti_sliceno(Roi, [])));
T1wName = 'T1W_brain_restore';
S_t1w = double(load_series(fullfile(Subject, T1wName), roi_nifti_sliceno(Roi, [])));

% Read Noise
fd = fopen(fullfile(Subject, 'GRE_restore_noise.txt'), 'r');
if fd ~= -1
    N_gre = double(fscanf(fd, '%f'));
    fclose(fd);
end

% Initialize output variables and start segmentation
N_iter = length(RoiLabelTable);
S_hypos_thr = zeros(size(S_roi), class(S_roi));
S_nontis = zeros(size(S_roi), class(S_roi));
S_ntis = zeros(size(S_roi), class(S_roi));
I_ntis_means = cell(N_iter, 1);
I_thr = cell(N_iter, 1);
CNR_oli = zeros(N_iter, 1);
for idx_iter = 1:N_iter
    Labs = RoiLabelTable{idx_iter};

    % Segmentation
    Ret = segment_us(S_gre, S_t1w, S_roi, S_ref, Labs, [], [], AdaptiveFlag, ...
                     'Out_name', ReportFile, 'P_thr', ThreshFactor);

    % Accumulate results
    S_hypos_thr = S_hypos_thr + Ret.S_hypos;
    S_nontis = S_nontis + Ret.S_nontis;
    S_ntis = S_ntis + Ret.S_ntis;
    I_ntis_means{idx_iter} = Ret.I_ntis_means';
    I_thr{idx_iter} = Ret.I_thr';
    CNR_oli(idx_iter) = abs(Ret.I_thr(1, 1)-Ret.I_ntis_means(1, 1))/N_gre;
end

% CC Filtering
[S_hypos, S_hypos_hypo, S_hypos_hyper, CC_app] = cc_filter(S_gre, S_t1w, ...
    N_gre, S_roi, logical(S_hypos_thr), S_ntis, [I_thr{:}]', ...
    CNR_thr, intvar_thr);

% Summary plot
H = figure; %create_ps_figure;
Col = [[0, 0, 0]; hsv(3)];
scatter(S_gre(logical(S_roi)), S_t1w(logical(S_roi)), 10, Col(1, :));
hold on;
I_ntis_means_list = [I_ntis_means{:}]';
scatter(I_ntis_means_list(:, 1), I_ntis_means_list(:, 2), 20, 'r');
SM_hypos = logical(S_hypos);
Mat = [S_gre(SM_hypos) S_t1w(SM_hypos)];
scatter(Mat(:, 1), Mat(:, 2), 10, Col(2, :));
if ~isempty(S_ref)
    SM_ref = logical(S_ref);
    scatter(S_gre(SM_ref & SM_hypos), S_t1w(SM_ref & SM_hypos), 10, Col(3, :));
    scatter(S_gre(SM_ref & ~SM_hypos), S_t1w(SM_ref & ~SM_hypos), 10, Col(4, :));
end
title(sprintf('Hypointense voxels: %d', sum(SM_hypos(:))));
save_ps_figure(ReportFile, H);

if SaveMaskFlag
    % Save masks
    save_series(RoiFile, fullfile(Subject, 'T2swHypo_mask'), ...
                S_hypos, roi_nifti_sliceno(Roi, []));
	save_series(RoiFile, fullfile(Subject, 'T2swHypoT1whypo_mask'), ...
                S_hypos_hypo, roi_nifti_sliceno(Roi, []));
	save_series(RoiFile, fullfile(Subject, 'T2swHypoT1whyper_mask'), ...
                S_hypos_hyper, roi_nifti_sliceno(Roi, []));
    S_tmp = S_hypos;
    S_tmp(logical(S_hypos_hypo)) = 0;
    S_tmp(logical(S_hypos_hyper)) = 0;
	save_series(RoiFile, fullfile(Subject, 'T2swHypoT1wiso_mask'), ...
                S_tmp, roi_nifti_sliceno(Roi, []));
    % Save non tissue mask
    save_series(RoiFile, fullfile(Subject, 'NonTis_mask'), ...
                S_nontis, roi_nifti_sliceno(Roi, []));
    % Save norm tissue mask
    save_series(RoiFile, fullfile(Subject, 'NormTis_mask'), ...
                S_ntis, roi_nifti_sliceno(Roi, []));
%     % Save delta delta R2 dash map
%     ddR2d = get_ddR2smap(S_gre, 15e-3, S_t1w, 102.96e-3, S_roi, I_ntis_means_list);
%     save_series(RoiFile, fullfile(Subject, 'ddR2d_map'), ...
%                 single(ddR2d), roi_nifti_sliceno(Roi, []));
end

% Validate
NII = load_series(RoiFile, 0);
F = NII.hdr.dime.pixdim(2:4);
if isempty(S_ref)
    S_ref = zeros(size(S_hypos));
end
Ret = validate_raw(S_hypos, S_ref, [], F);
Ret.thr = validate_raw(S_hypos_thr, S_ref, [], F);

% Validate hand-edited masks
% try
%     FeName = 'FE_roi_mask_edited';
%     S_ref_orig = load_series(fullfile(Subject, FeName), roi_nifti_sliceno(Roi, []));
%     S_ref_orig(S_ref_orig == 50) = 11;
%     S_ref_orig(S_ref_orig == 51) = 12;
%     S_ref_orig(S_ref_orig == 52) = 13;
%     S_ref_orig(S_ref_orig == 55) = 14;
%     S_ref_edit = zeros(size(S_ref_orig), class(S_ref_orig));
%     N_iter = length(RoiLabelTable);
%     for idx_iter = 1:N_iter
%         Labs = RoiLabelTable{idx_iter};
%         for lab = Labs
%             S_ref_edit(S_ref_orig == lab) = lab;
%         end
%     end
% catch
%     S_ref_edit = S_hypos;
% end
% Ret.edit = validate_raw(S_ref_edit, S_ref, [], F);

% Add input parameters to output
Ret.Input.Subject = Subject;
Ret.Input.RoiLabelTable = RoiLabelTable;
Ret.Input.ReportName = ReportName;
Ret.Input.ThreshFactor = ThreshFactor;
Ret.Input.AdaptiveFlag = AdaptiveFlag;
Ret.Input.N_gre = N_gre;
Ret.Input.CNR_thr = CNR_thr;
Ret.Input.intvar_thr = intvar_thr;

% Add accumulated results to output
Ret.I_thr = I_thr;
Ret.I_ntis_means = I_ntis_means;
Ret.CNR_oli = CNR_oli;
Ret.CC_app = CC_app;

