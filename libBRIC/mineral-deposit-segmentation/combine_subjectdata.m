% First include our libraries
addpath('../../libBRIC/misc-matlab/');

% Include external libraries
addpath('NIFTI/');

close all; clear all; clc;

anot_file = '/media/LP3TBdisk/Andreas_PhD/mineral/annot.csv';
base_path = '/media/LP3TBdisk/Andreas_PhD/mineral-deposit-segmentation-pipeline/BRICpipe';
base_name = 'subjects_asps_318_checked_both_t2w';
sub_file = fullfile(base_path, [base_name '.csv']);
mat_name = [base_name '_all.mat']; % matlab file for saving results
ps_name = [base_name '_all'];
t2swhypo_name = 'T2swHypo_mask'; % name of T2*w hypointensity mask subject directory
t2swhypot1whypo_name = 'T2swHypoT1whypo_mask';
t2swhypot1whyper_name = 'T2swHypoT1whyper_mask';
ntis_name = 'NormTis_mask';
roi_name = 'RO_weinew_mask'; % name of ROI mask in subject directory
t2sw_name = 'GRE_brain_restore'; % name of GRE volume
t1w_name = 'T1W_brain_restore';
t2w_name = 'T2W_brain_restore';
r2s_name = 'R2s';
r2_name = 'R2';
r2d_name = 'R2d';
pha_name = 'GRE_phu';
pha_valid_name = 'GRE_sos_bin';
icv_name = 'GRE_brain_mask'; % name of ICV mask

% Read subject file
[Subjects, N_total] = get_subject_dirs(sub_file);

% Read annot file
fd = fopen(anot_file);
SubjAnot = textscan(fd, '%s%s%s%s%s%s%*[^\n]', ...
    'delimiter',',',...
    'treatAsEmpty',{'NA','na'}, ...
    'commentStyle', '#');
fclose(fd);

% Init result structure
Ret = struct;
ret_idx = 1;
for i=1:N_total
    path = char(Subjects(i));
    Ret(ret_idx).Path = path;
    fprintf('Reading subject data for subject %s...', path);
    
    S_feat = load_series(fullfile(path, t2swhypo_name), []);
    N_fe = sum(logical(S_feat(:)));
    fprintf('%d voxels...', N_fe);
    
    % Find corresponding annotation (chi number that matches path)
    f = @(in) (~isempty(strfind(path, in)));
    M_subj = cellfun(f, SubjAnot{1});
    
    % Skip empty masks
    if N_fe > 0 && sum(M_subj) == 1
        % Get connected component statistics
        [Ret(ret_idx).CC] = ...
            conncomp_stats(path, t2swhypo_name, 0, roi_name, 1, 1, ...
                           t2sw_name, t2swhypot1whypo_name, t2swhypot1whyper_name);

        % Labels of ROI mask
        S_roi = load_series(fullfile(path, roi_name), []);
        Lab_roi = unique(S_roi(:))';
        Ret(ret_idx).Lab_roi = Lab_roi(2:end); % creates list of roi labels
        
        % Labels of ID mask
        S_feat = load_series(fullfile(path, t2swhypo_name), []);
        Ret(ret_idx).Lab = S_feat(logical(S_feat))'; % creates a line vector
        
        % ICV
        SM_brain = logical(load_series(fullfile(path, icv_name), []));
        Ret(ret_idx).ICV = str2double(SubjAnot{5}{M_subj});
        if Ret(ret_idx).ICV <= 0 || isnan(Ret(ret_idx).ICV)
            NII_brain = load_series(fullfile(path, icv_name), 0);
            Ret(ret_idx).ICV = get_volume(SM_brain, ...
                                          NII_brain.hdr.dime.pixdim(2:4));
        end
        
        % Age
        Ret(ret_idx).Age = str2double(SubjAnot{2}{M_subj});
        
        % Gender
        Ret(ret_idx).Male = strcmp(SubjAnot{6}{M_subj}, 'm');
        Ret(ret_idx).Female = strcmp(SubjAnot{6}{M_subj}, 'w');

        % T2sw Intensity info
        S_t2sw = single(load_series(fullfile(path, t2sw_name), []));
        S_ntis = single(load_series(fullfile(path, ntis_name), []));
        N_lab = length(Ret(ret_idx).Lab_roi);
        for lab_idx = 1:N_lab
            SM_ntis = S_ntis == Ret(ret_idx).Lab_roi(lab_idx);
            SM_feat = S_feat == Ret(ret_idx).Lab_roi(lab_idx);
            SM_valid = S_t2sw > 0 & SM_brain;
            Ret(ret_idx).T2sw(lab_idx) = get_tis_ints(S_t2sw, SM_ntis, SM_feat, SM_valid);
        end
        Ret(ret_idx).T2swValid = 1;
        
        % T1w Intensity info
        S_t1w = single(load_series(fullfile(path, t1w_name), []));
        N_lab = length(Ret(ret_idx).Lab_roi);
        for lab_idx = 1:N_lab
            SM_ntis = S_ntis == Ret(ret_idx).Lab_roi(lab_idx);
            SM_feat = S_feat == Ret(ret_idx).Lab_roi(lab_idx);
            SM_valid = S_t1w > 0 & SM_brain;
            Ret(ret_idx).T1w(lab_idx) = get_tis_ints(S_t1w, SM_ntis, SM_feat, SM_valid);
        end
        Ret(ret_idx).T1wValid = 1;
        
        % T2w Intensity info
        S_t2w = single(load_series(fullfile(path, t2w_name), []));
        N_lab = length(Ret(ret_idx).Lab_roi);
        for lab_idx = 1:N_lab
            SM_ntis = S_ntis == Ret(ret_idx).Lab_roi(lab_idx);
            SM_feat = S_feat == Ret(ret_idx).Lab_roi(lab_idx);
            SM_valid = S_t2w > 0 & SM_brain;
            Ret(ret_idx).T2w(lab_idx) = get_tis_ints(S_t2w, SM_ntis, SM_feat, SM_valid);
        end
        Ret(ret_idx).T2wValid = 1;
        
        % R2s Intensity info
        S_r2s = single(load_series(fullfile(path, r2s_name), []));
        S_r2s = S_r2s(:, :, :, 1);
        N_lab = length(Ret(ret_idx).Lab_roi);
        for lab_idx = 1:N_lab
            SM_ntis = S_ntis == Ret(ret_idx).Lab_roi(lab_idx);
            SM_feat = S_feat == Ret(ret_idx).Lab_roi(lab_idx);
            SM_valid = S_r2s > 0 & SM_brain;
            Ret(ret_idx).R2s(lab_idx) = get_tis_ints(S_r2s, SM_ntis, SM_feat, SM_valid);
        end
        Ret(ret_idx).R2sValid = 1;
        
        % R2 Intensity info
        S_r2 = single(load_series(fullfile(path, r2_name), []));
        N_lab = length(Ret(ret_idx).Lab_roi);
        for lab_idx = 1:N_lab
            SM_ntis = S_ntis == Ret(ret_idx).Lab_roi(lab_idx);
            SM_feat = S_feat == Ret(ret_idx).Lab_roi(lab_idx);
            SM_valid = S_r2 > 0 & SM_brain;
            Ret(ret_idx).R2(lab_idx) = get_tis_ints(S_r2, SM_ntis, SM_feat, SM_valid);
        end
        Ret(ret_idx).R2Valid = 1;
        
                
        % R2d Intensity info
        S_r2d = single(load_series(fullfile(path, r2d_name), []));
        N_lab = length(Ret(ret_idx).Lab_roi);
        for lab_idx = 1:N_lab
            SM_ntis = S_ntis == Ret(ret_idx).Lab_roi(lab_idx);
            SM_feat = S_feat == Ret(ret_idx).Lab_roi(lab_idx);
            SM_valid = S_r2d > 0 & SM_brain;
            Ret(ret_idx).R2d(lab_idx) = get_tis_ints(S_r2d, SM_ntis, SM_feat, SM_valid);
        end
        Ret(ret_idx).R2dValid = 1;
        
        try
            % Pha Intensity info
            S_pha = single(load_series(fullfile(path, pha_name), []));
            S_pha = S_pha(:, :, :, 3);
            SM_pha_valid = logical(load_series(fullfile(path, pha_valid_name), []));
            N_lab = length(Ret(ret_idx).Lab_roi);
            for lab_idx = 1:N_lab
                SM_ntis = S_ntis == Ret(ret_idx).Lab_roi(lab_idx);
                SM_feat = S_feat == Ret(ret_idx).Lab_roi(lab_idx);
                SM_valid = SM_pha_valid & SM_brain;
                Ret(ret_idx).Pha(lab_idx) = get_tis_ints(S_pha, SM_ntis, SM_feat, SM_valid);
            end
            Ret(ret_idx).PhaValid = 1;
        catch
            Ret(ret_idx).Pha = Ret(ret_idx).R2;
            Ret(ret_idx).PhaValid = 0;
        end
            
        ret_idx = ret_idx + 1;
    end
    
    fprintf('done.\n');
end

save(mat_name);
% Uncomment from here upwards if run once

load(mat_name);
% M = [Ret.Age] < 65;
% Ret = Ret(M);
% ps_name = [base_name '_iso_lt65'];

% Summary stats
N_subject = length(Ret);
N_CC_roi = get_cc_avg_stats({[Ret.CC]}, Ret(1).Lab_roi, [], '', '');
N_CC = zeros(size(Ret(1).Lab_roi));
N_roi = length(Ret(1).Lab_roi);
for idx = 1:N_roi
    N_CC(idx) = get_cc_avg_stats({[Ret.CC]}, Ret(1).Lab_roi(idx), [], '', '');
end
fprintf('- Summary -------------------------------------------------------\n');
fprintf('#Subjects: %d\n', N_subject);
fprintf('#CC roi: %d\n', N_CC_roi);
fprintf('#CC:'); [Ret(1).Lab_roi; N_CC], sum(N_CC)
fprintf('=================================================================\n');

save_ps_figure(ps_name, []);
% Intensities
[H] = print_intvsastructure(Ret, 'T2sw', 3, 0.75, 1);
save_ps_figure(ps_name, H);
[H] = print_intvsastructure(Ret, 'T2w', 3, 0.75, 1);
save_ps_figure(ps_name, H);
[H] = print_intvsastructure(Ret, 'T1w', 3, 0.75, 1);
save_ps_figure(ps_name, H);
[H] = print_intvsastructure(Ret, 'R2s', 3, 0.75);
save_ps_figure(ps_name, H);
[H] = print_intvsastructure(Ret, 'R2', 3, 0.75);
save_ps_figure(ps_name, H);
[H] = print_intvsastructure(Ret, 'R2d', 3, 0.75);
save_ps_figure(ps_name, H);
M = [Ret.PhaValid] > 0;
[H] = print_intvsastructure(Ret(M), 'Pha', 3, 0.75);
save_ps_figure(ps_name, H);
% Appearance
[H] = print_appvsastructure(Ret, 'T2sw', 0.75);
save_ps_figure(ps_name, H);
[H] = print_appvsastructure(Ret, 'T2w', 0.75);
save_ps_figure(ps_name, H);
[H] = print_appvsastructure(Ret, 'T1w', 0.75);
save_ps_figure(ps_name, H);
% [H] = print_appvsastructure(Ret, 'R2s');
% save_ps_figure(ps_name, H);
% [H] = print_appvsastructure(Ret, 'R2');
% save_ps_figure(ps_name, H);
M = [Ret.PhaValid] > 0;
[H] = print_appvsastructure(Ret(M), 'Pha', 0.75);
save_ps_figure(ps_name, H);

% Per subject count & volume
Loc = Ret(1).Lab_roi;
LocMat = zeros(N_subject, length(Loc)+1);
VolMat = zeros(N_subject, length(Loc)+1);
for idx = 1:length(Loc)
    [LocMat(:, idx), Val]  = get_cc_avg_stats({Ret.CC}, Loc(idx), 'vol', '', 'sum');
    VolMat(:, idx) = cell2mat(Val);
end
LocMat(:, length(Loc)+1) = sum(LocMat(:, 1:4), 2);
VolMat(:, 1:length(Loc)) = VolMat(:, 1:length(Loc)) ./ repmat([Ret.ICV]', 1, length(Loc));
VolMat(:, length(Loc)+1) = sum(VolMat(:, 1:4), 2);
H = figure; %create_ps_figure;
subplot(2,1,1);
hl = boxplot(LocMat+eps, 'plotstyle', 'compact'); %'notch', 'on');
for ih=1:size(hl, 1)
    set(hl(ih, :), 'Color', 'k');
end
set(gca, 'XTick', 1:size(LocMat, 2));
set(gca, 'XTickLabel', {'C','P', 'G', 'I', 'Left+Right'}');
xlabel('\bf Anatomical structure(s)');
ylabel('\bf BGID count per subject');
fprintf('-- Relative BGID count per subject --\n');
quantile(LocMat, [.25 .5 .75])
fprintf('-------------------------------------\n');
subplot(2,1,2);
hl = boxplot(VolMat.*1e6+eps, 'plotstyle', 'compact'); %'notch', 'on');'notch', 'on');
for ih=1:size(hl, 1)
    set(hl(ih, :), 'Color', 'k');
end
set(gca, 'XTick', 1:size(LocMat, 2));
set(gca, 'XTickLabel', {'C','P', 'G', 'I', 'Left+Right'}');
xlabel('\bf Anatomical structure(s)');
ylabel(sprintf('\\bfBGID load per subject\n\\bfin ppm of ICV'));
fprintf('-- Relative BGID volume per subject in ppm of ICV --\n');
quantile(VolMat.*1e6, [.25 .5 .75])
fprintf('---------------------------------------------------------\n');
save_ps_figure(ps_name, H);

% Total volume
H = figure; %create_ps_figure;
Tmp = sum(VolMat(:, 1:length(Loc)), 1).*1000;
bar(Tmp);
xlabel('\bf Location');
ylabel('\bf Total relative BGID volume in permille of ICV');
%save_ps_figure(OutFile, H);
fprintf('-- Total Relative BGID volume in permille of ICV --\n');
Tmp
fprintf('---------------------------------------------------\n');
save_ps_figure(ps_name, H);

% Stats of CC within an area (vol, max. area, compactness, rel. aniso)
Loc = cell(3, 1);
Loc{1} = Ret(1).Lab_roi([1 2 4]);
Loc{2} = Ret(1).Lab_roi([3]);
Loc{3} = Ret(1).Lab_roi;
Q = [.25 .5 .75];
Field = cell(4, 1);
Field{1} = 'vol';
Field{2} = 'ma';
Field{3} = 'com';
Field{4} = 'ra';
Field{5} = 'phypo';
Field{6} = 'phyper';
CC_Type = cell(3, 1);
CC_Type{1} = '';
CC_Type{2} = 'inter';
CC_Type{3} = 'intra';
for idx = 1:length(CC_Type)
    [ResMat, N_res] = get_cc_avg_stats_all({Ret.CC}, Loc, Q, Field, CC_Type{idx});
    fprintf(['-- CC ' CC_Type{idx} ' stats (vol, ma, com, ra, phypo, phyper) --\n']);
    N_res'
    fprintf('C, P, I'); reshape(ResMat(:, 1, :), 6, 3)
    fprintf('------------------------------------------------------------\n');
    fprintf('G'); reshape(ResMat(:, 2, :), 6, 3)
    fprintf('------------------------------------------------------------\n');
    fprintf('All'); reshape(ResMat(:, end, :), 6, 3)
    fprintf('============================================================\n');
end
