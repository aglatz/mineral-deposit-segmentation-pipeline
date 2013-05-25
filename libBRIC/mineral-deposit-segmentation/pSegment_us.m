% This script gets called from segment_us_mp(), which distributes the
% input data accross several CPU cores using pMatlab (see
% http://www.ll.mit.edu/mission/isr/pmatlab/ for further info).
%
% This script expects three files in each subject directory:
% - RO_mask: the masks with the regions of interests
% - GRE_brain_restore: the GRE volume which was brain extracted
%                      (not strictly necessary) and bias-field corrected.
% - T1W_brain_restore: the T1W volume which was brain extracted
%                      (not strictly necessary) and bias-field corrected
% It generates the following files for each subject:
% - NonTis_mask: A mask selecting what appears to be not normal tissue
% - NormTis_mask: A mask selecting the normal-appearing tissue
% - T2swHypo_mask: A mask selecting T2sw hypointensities (subset of
%                  NonTis_mask)
%

% First include our libraries
addpath('../misc-matlab/');

% Include external libraries - create symbolic links if it fails here!
addpath('NIFTI/');
addpath('LIBRA/');

close all; % No 'clear all' otherwise we loose our input variables!

% Load input vars in mp mode
InputFileVar = 'InputFile';
if exist(InputFileVar, 'var') > 0
    PMatDir = tempdir;
    InputFile = [PMatDir '/' InputFile];
    load(InputFile);
end

% Default input arguments
if exist('ReportName', 'var') <= 0 || isempty(ReportName)
    ReportName = 'class';
end
if exist('InterpFactor', 'var') <= 0 || isempty(InterpFactor)
    InterpFactor = 1;
end
if exist('ThreshFactor', 'var') <= 0 || isempty(ThreshFactor)
    ThreshFactor = [0.9579 0.4314];
end

% Read subject file - last column contains the paths to the subject dirs
[ndata, text, raw] = xlsread(SubjectFile); clear ndata text;
Subjects = raw(:, size(raw, 2));
N_total = length(Subjects);

if exist(InputFileVar, 'var') > 0
    % pMatlab setup - Init shared matrix
    Wmap=map([Np 1], {}, 0:Np-1);
    % Wmap=1; % Uncomment for sequential processing
    Over = zeros(N_total, length(RoiLabelTable)*2, Wmap);
    % Get portion of matrix for this process
    OverLoc = local(Over);
    IdxLoc = global_ind(Over, 1);
else
    OverLoc = zeros(N_total, length(RoiLabelTable)*2, 1);
    IdxLoc = 1:N_total;
end
for idx_j = 1:size(OverLoc, 1);
	idx_sub = IdxLoc(idx_j);
	Subject = Subjects{idx_sub};
    fprintf('Subject: %s ...\n', Subject);
    
    ReportFile = [Subject '/' ReportName];
    save_ps_figure(ReportFile, []); % Deletes previous file

    % Read ROI
    RoiName = 'RO_mask';
    RoiFile = [Subject '/' RoiName];
    S_roi = load_series(RoiFile, []);
    Roi = roi_init(S_roi);
    S_roi = load_series_interp(RoiFile, roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);
    
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
    S_nontis_all = zeros(size(S_roi), class(S_roi));
    S_ntis_all = zeros(size(S_roi), class(S_roi));
    I_ntis_means = cell(N_iter, 1);
    for idx_iter = 1:N_iter
        Labs = RoiLabelTable{idx_iter};
        N_labs = length(Labs);
        
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
        [S_out, S_nontis, S_ntis, I_ntis_means_iter, I_gre_thr] = segment_us(...
                                S_gre, S_t1w, S_roi, S_ref, Labs, [], [], ...
                                'Out_name', ReportFile, 'P_thr', ThreshFactor);
                   
        % Store results
        S_out_all = S_out_all + S_out;
        S_nontis_all = S_nontis_all + S_nontis;
        S_ntis_all = S_ntis_all + S_ntis;
        I_ntis_means{idx_iter} = I_ntis_means_iter;
        Idx_over = (idx_iter-1)*N_iter + [1:2];
        OverLoc(idx_j, Idx_over(1)) = I_gre_thr;
        if ~isempty(S_ref)
            SM_tmp1 = false(size(S_out));
            for lab = Labs
                SM_tmp1 = SM_tmp1 | S_out == lab;
            end
            SM_tmp2 = false(size(S_ref));
            for lab = Labs
                SM_tmp2 = SM_tmp2 | S_ref == lab;
            end
            OverLoc(idx_j, Idx_over(2)) = quantile(S_gre(SM_tmp1 & SM_tmp2), .95); 
        end
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
    
    % Save hypomask
	T2swHypoMaskName = 'T2swHypo_mask';
    S_out_all = interp_series([Subject '/' RoiName], S_out_all, ...
                    roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);
    [SM_oli] = morph_filter([Subject '/' T2swName], ...
                        roi_nifti_sliceno(Roi, []), logical(S_out_all));
    S_out_all(~SM_oli) = 0;
    save_series([Subject '/' RoiName], [Subject '/' T2swHypoMaskName], ...
                S_out_all, roi_nifti_sliceno(Roi, []));
    % Save non tissue mask
    save_series_interp([Subject '/' RoiName], [Subject '/NonTis_mask'], ...
                        S_nontis_all, roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);
    % Save norm tissue mask
    save_series_interp([Subject '/' RoiName], [Subject '/NormTis_mask'], ...
                        S_ntis_all, roi_nifti_sliceno(Roi, []), 'nearest', InterpFactor);
    
	% Validate
    if ~isempty(S_ref)
        Ret = validate(Subject, T2swHypoMaskName, FeName, [], RoiName, 1);
        save([Subject '/Ret.mat'], 'Ret');
    else
        % Fake validation just to get the volume
        Ret = validate(Subject, T2swHypoMaskName, T2swHypoMaskName, [], RoiName, 1);
        save([Subject '/Ret.mat'], 'Ret');
    end
end

if exist(InputFileVar, 'var') > 0
    % pMatlab exit - Aggregate shared matrix
    Over = put_local(Over, OverLoc);
    Over = agg(Over);
else
    Over = OverLoc;
end

% Collect result data in parent process (Pid=0)
if ~Pid && ~isempty(S_ref)
    % Combine all returned structures to a array of structures
    RetAll = Ret(end); % just inititializes the structure array
    for i = 1:N_total
        RetPath = [char(Subjects(i)) '/Ret.mat'];
        load(RetPath);
        RetAll(i) = Ret(end); % Last one is best one
    end
end
