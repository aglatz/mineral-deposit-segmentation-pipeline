% First include our libraries
addpath('../../libBRIC/misc-matlab/');

% Include external libraries
addpath('NIFTI/');

close all; clear all; clc;

sub_file = 'subjects_80_all.csv'; % csv file
mat_name = 'subjects_80_all.mat'; % matlab file for saving results
fe_name = 'FE_roi_mask'; % name of fe mask subject directory
roi_name = 'RO_mask'; % name of ROI mask in subject directory
gre_name = 'GRE_brain_restore'; % name of GRE volume
icv_name = 'GRE_brain_mask'; % name of ICV mask

% Read subject file
fd = fopen(sub_file);
if fd < 0
    error('create_subjectdata:fileopenerror', 'Could not open file!');
end
Subjects = textscan(fd, '%s%s%s%s%s%s%s%s%s%*[^\n]', ...
                       'delimiter',',',...
                       'treatAsEmpty',{'NA','na'}, ...
                       'commentStyle', '#');
fclose(fd);
Subjects = Subjects{9}; % contains path to subject directories
N_total = length(Subjects);

% Init result structure
Ret = struct;
ret_idx = 1;
for i=1:N_total
    path = char(Subjects(i));
    Ret(ret_idx).Path = path;
    fprintf('Reading subject data for subject %s...', path);
    
    S_fe = load_series(fullfile(path, fe_name), []);
    N_fe = sum(logical(S_fe(:)));
    fprintf('%d voxels...', N_fe);
    
    % Skip empty masks
    if N_fe > 0
        % Get connected component statistics
        [Ret(ret_idx).CC] = ...
            conncomp_stats(path, fe_name, 0, roi_name, 1, 1, gre_name);

        % Labels of ROI mask
        S_roi = load_series(fullfile(path, roi_name), []);
        Lab_roi = unique(S_roi(:))';
        Ret(ret_idx).Lab_roi = Lab_roi(2:end); % creates list of roi labels
        
        % Labels of ID mask
        S_fe = load_series(fullfile(path, fe_name), []);
        Ret(ret_idx).Lab = S_fe(logical(S_fe))'; % creates a line vector
        
        % ICV
        SM_brain = logical(load_series(fullfile(path, icv_name), []));
        NII_brain = load_series(fullfile(path, icv_name), 0);
        Ret(ret_idx).ICV = get_volume(SM_brain, ...
                                NII_brain.hdr.dime.pixdim(2:4));
                            
        ret_idx = ret_idx + 1;
    end
    
    fprintf('done.\n');
end

save(mat_name);
% Uncomment from here upwards if run once

load(mat_name);

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

% Per subject count & volume
Loc = Ret(1).Lab_roi;
LocMat = zeros(N_subject, length(Loc)+3);
VolMat = zeros(N_subject, length(Loc)+3);
for idx = 1:length(Loc)
    [LocMat(:, idx), Val]  = get_cc_avg_stats({Ret.CC}, Loc(idx), 'vol', '', 'sum');
    VolMat(:, idx) = cell2mat(Val);
end
LocMat(:, length(Loc)+1) = sum(LocMat(:, 1:4), 2);
LocMat(:, length(Loc)+2) = sum(LocMat(:, 5:length(Loc)), 2);
LocMat(:, length(Loc)+3) = sum(LocMat(:, 1:length(Loc)), 2);
VolMat(:, 1:length(Loc)) = VolMat(:, 1:length(Loc)) ./ repmat([Ret.ICV]', 1, length(Loc));
VolMat(:, length(Loc)+1) = sum(VolMat(:, 1:4), 2);
VolMat(:, length(Loc)+2) = sum(VolMat(:, 5:length(Loc)), 2);
VolMat(:, length(Loc)+3) = sum(VolMat(:, 1:length(Loc)), 2);
H = figure; %create_ps_figure;
subplot(2,1,1);
hl = boxplot(LocMat+eps, 'plotstyle', 'compact'); %'notch', 'on');
for ih=1:size(hl, 1)
    set(hl(ih, :), 'Color', 'k');
end
set(gca, 'XTick', 1:size(LocMat, 2));
set(gca, 'XTickLabel', {'CL','PL', 'GL', 'IL', ...
                        'CR', 'PR', 'GR', 'IR', 'Left', 'Right', 'Left+Right'}');
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
set(gca, 'XTickLabel', {'CL','PL', 'GL', 'IL', ...
                        'CR', 'PR', 'GR', 'IR', 'Left', 'Right', 'Left+Right'}');
xlabel('\bf Anatomical structure(s)');
ylabel(sprintf('\\bfBGID load per subject\n\\bfin ppm of ICV'));
fprintf('-- Relative BGID volume per subject in ppm of ICV --\n');
quantile(VolMat.*1e6, [.25 .5 .75])
fprintf('---------------------------------------------------------\n');
%save_ps_figure(OutFile, H);

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

% Stats of CC within an area (vol, max. area, compactness, rel. aniso)
Loc = cell(3, 1);
Loc{1} = Ret(1).Lab_roi([1 2 4 5 6 8]);
Loc{2} = Ret(1).Lab_roi([3 7]);
Loc{3} = Ret(1).Lab_roi;
Q = [.25 .5 .75];
Field = cell(4, 1);
Field{1} = 'vol';
Field{2} = 'ma';
Field{3} = 'com';
Field{4} = 'ra';
CC_Type = cell(3, 1);
CC_Type{1} = '';
CC_Type{2} = 'inter';
CC_Type{3} = 'intra';
for idx = 1:length(CC_Type)
    [ResMat, N_res] = get_cc_avg_stats_all({Ret.CC}, Loc, Q, Field, CC_Type{idx});
    fprintf(['-- CC ' CC_Type{idx} ' stats (vol, max. area, compactness, rel. anisotropy --\n']);
    N_res'
    fprintf('C, P, I'); reshape(ResMat(:, 1, :), 4, 3)
    fprintf('------------------------------------------------------------\n');
    fprintf('G'); reshape(ResMat(:, 2, :), 4, 3)
    fprintf('------------------------------------------------------------\n');
    fprintf('All'); reshape(ResMat(:, end, :), 4, 3)
    fprintf('============================================================\n');
end
