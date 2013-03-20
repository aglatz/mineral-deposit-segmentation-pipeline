% This script gets called from segment_us_mp(), which distributes the
% input data accross several CPU cores using pMatlab (see
% http://www.ll.mit.edu/mission/isr/pmatlab/ for further info). To run
% this script standalone (un)comment lines below.
% This script expects three files in each subject directory:
% - RO_mask: the masks with the regions of interests
% - GRE_brain_restore: the GRE volume which was brain extracted
%                      (not strictly necessary) and bias-field corrected.
% - T1W_brain_restore: the T1W volume which was brain extracted
%                      (not strictly necessary) and bias-field corrected
% If a reference mask is available for a subject then the mask should
% have the file name 'FE_roi_mask'.
% INPUTS (as global variables): SubjectFile - A csv file in the format
%                       used for preprocessing. The last column is read
%                       and should contain the subject directory paths.
% RETURNS (as global variables): Ret - Cell array containing the results
%                                      from validate() for every subject.
%                                Over - Contains info about GRE cutoff
%                       thresholds calculated with the generated and
%                       reference masks, if given. See last two return
%                       values of segment_us().
%                                Subjects - Cell array with subjects paths.
%

addpath('NIFTI/');
addpath('LIBRA/');
addpath('libBRIC/');

close all;

if isempty(SubjectFile)
	error('pSegment:Inputargs', 'Insufficient # of input arguments');
end
OutFile = 'class';

% Read subject file
fd = fopen(SubjectFile);
Subjects = textscan(fd, '%s%s%s%s%s%s%s%s%s%*[^\n]', ...
                       'delimiter',',',...
                       'treatAsEmpty',{'NA','na'}, ...
                       'commentStyle', '#');
fclose(fd);
Subjects = Subjects{9};
% load('subjects_65.mat');
N_total = length(Subjects);

% Init shared matrix
Wmap=map([Np 1], {}, 0:Np-1);
% Wmap=1; % Uncomment for sequential processing
Over = zeros(N_total, 12, Wmap);

% Get portion of matrix for this process
OverLoc = local(Over);
IdxLoc = global_ind(Over, 1);
for idx_j = 1:size(OverLoc, 1);
	idx_sub = IdxLoc(idx_j);
	Subject = Subjects{idx_sub};
    fprintf('Subject: %s ...\n', Subject);
    
    Out_name = [Subject '/' OutFile];
    save_ps_figure(Out_name, []); % Delete previous file

    % Read data
    Name_roi_red = 'RO_mask';
    Name_roi = [Subject '/' Name_roi_red];
    S_roi = load_series(Name_roi, []);
    Roi = roi_init(S_roi);
    S_roi = load_series(Name_roi, roi_nifti_sliceno(Roi, []));
    try
        Fe_name = 'FE_roi_mask';
        S_ref = load_series([Subject '/' Fe_name], roi_nifti_sliceno(Roi, []));
    catch
        S_ref = [];
    end
    GRE_name = 'GRE_brain_restore';
    S_gre = double(load_series([Subject '/' GRE_name], roi_nifti_sliceno(Roi, [])));
    T1W_name = 'T1W_brain_restore';
    S_t1w = double(load_series([Subject '/' T1W_name], roi_nifti_sliceno(Roi, [])));
    
% 	SM_wm = logical(load_series([Subject '/WM_mask'], roi_nifti_sliceno(Roi, [])));
%     [~, I_gre_wm_mu] = volume_stats(S_gre, SM_wm, [], [], []);
%     TE = 15e-3;
%     I_gre_thr = I_gre_wm_mu * exp(-15e-3*12.35);
%     I_gre_thr = [];
    
    % Segment
    Lab = unique(S_roi)';
    Lab = Lab(2:end); % Exclude background
    Lab = Lab([3 1 2 4 7 5 6 8]);
    N_lab = length(Lab);
    S_out_all = zeros(size(S_roi), class(S_roi));
    S_nontis_all = zeros(size(S_roi), class(S_roi));
    S_ntis_all = zeros(size(S_roi), class(S_roi));
    I_means_est = zeros(N_lab, 2);
    for idx_step = 1:2
        idx_off = 4 * (idx_step-1);
        Idx_lab = [1 2 3 4] + idx_off;
        
        % Heuristic to estimate the normal-appearing tissue means - works
        % for GRE and T1W of LBC1936 protocol
        I_ntis_means_all = NaN(length(Idx_lab), 2);
        for idx_lab = 1:length(Idx_lab)
            SM_tmp = S_roi==Lab(Idx_lab(idx_lab));
            Mat = [S_gre(SM_tmp) S_t1w(SM_tmp)];
            [~, I_ntis_means_all(idx_lab, :)] = pcomp_find(Mat);
        end
        Tmp = mcdregres(I_ntis_means_all(:, 2), I_ntis_means_all(:, 1), 'plots', 0);
        Est_param = [Tmp.slope Tmp.int];
        I_ntis_means_all_est = I_ntis_means_all;
        I_ntis_means_all_est(:, 1) = polyval(Est_param, I_ntis_means_all(:, 2));
        for idx_lab = 1:length(Idx_lab)
            idx_cur = Idx_lab(idx_lab);
            I_means_est(idx_cur, :) = I_ntis_means_all_est(idx_lab, :);
        end
        
        % Segmentation
        [S_out, S_nontis, S_ntis, I_ntis_means_all(Idx_lab, :), I_gre_thr, I_gre_thr_ref] = ...
            segment_us(S_gre, S_t1w, S_roi, S_ref, Lab(Idx_lab), I_means_est(Idx_lab, :), ...
                       [], 'Out_name', Out_name, 'P_thr', [0.9579 0.4314]);
        S_out_all = S_out_all + S_out;
        S_nontis_all = S_nontis_all + S_nontis;
        S_ntis_all = S_ntis_all + S_ntis;
        OverLoc(idx_j, (idx_step-1)*6+1) = I_gre_thr;
        OverLoc(idx_j, (idx_step-1)*6+2:(idx_step-1)*6+6) = I_gre_thr_ref;
    end
    
    % Summary plot
	H = create_ps_figure; %figure;
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
	Out_name = 'Oli_us_mask';
    save_series_roi(Subject, Name_roi_red, Out_name, S_out_all, roi_nifti_sliceno(Roi, []));
    save_series_roi(Subject, Name_roi_red, 'NonTis_mask', S_nontis_all, roi_nifti_sliceno(Roi, []));
    save_series_roi(Subject, Name_roi_red, 'NormTis_mask', S_ntis_all, roi_nifti_sliceno(Roi, []));
    
    preproc([], 1, Subject, T1W_name, GRE_name, ...
            'WM_mask', 'GM_mask', 'CS_mask', Out_name, Name_roi_red, ...
            'AR_mask', []);
    
    if ~isempty(S_ref)
        Ret = validate(Subject, Out_name, Fe_name, [], Name_roi_red, 1);
        save([Subject '/Ret.mat'], 'Ret');
    else
        Ret = validate(Subject, Out_name, Out_name, [], Name_roi_red, 1);
        save([Subject '/Ret.mat'], 'Ret');
    end
end

% Aggregate shared matrix
Over = put_local(Over, OverLoc);
Over = agg(Over);

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
