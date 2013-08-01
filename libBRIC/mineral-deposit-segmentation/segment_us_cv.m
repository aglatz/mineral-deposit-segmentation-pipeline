% This script first estimates the thresholds that lead to the best
% (or at least better) Jaccard indices, which are then used to
% refine the subject thresholds in a 10-fold cross-validation.
%

close all; clear all

SubjectFile = 'subjects_98';
Pfx = 'ad_100713';
AdaptiveFlag = true;
RoiLabelTable = {[13, 11, 12, 14]};
N_cpus = 6;

% Get optimal values
Ret = segment_us_mp([SubjectFile '.xls'], ...
                    RoiLabelTable, N_cpus, 'ThreshFactor', [1 0], ...
                    'FuncName', 'segment_us_refine', ...
                    'AdaptiveFlag', AdaptiveFlag);
save([SubjectFile '_' Pfx '.mat']);
load([SubjectFile '_' Pfx '.mat']);

% CV
[~, ~, raw] = xlsread([SubjectFile '.xls']);
Subjects = raw(:, size(raw, 2));
N_rep = 1;
N_total = length(Subjects);
for idx_rep = 1:N_rep
    cvp = cvpartition(N_total, 'k', 10);
    out = struct;
    for idx_testset = 1:cvp.NumTestSets
        % Train
        idx_train = cvp.training(idx_testset);
        Alpha = [Ret(idx_train).alpha_opt];
        M = ~isnan(Alpha); % Ignore non optimized subjects
        out(idx_testset).P = [median(Alpha(M)) 0];

        % Test
        idx_test = cvp.test(idx_testset);
        Subjects_test = Subjects(idx_test);
        SubjectFile_traintest = [SubjectFile '_' Pfx '_' num2str(idx_testset)];
        save_xls(SubjectFile_traintest, Subjects_test);
        [out(idx_testset).Ret, out(idx_testset).Subjects] = ...
                        segment_us_mp([SubjectFile_traintest '.xls'], ...
                                      RoiLabelTable, N_cpus, ...
                                      'ThreshFactor', out(idx_testset).P, ...
                                      'ReportName', 'class', ...
                                      'AdaptiveFlag', AdaptiveFlag);
    end
    save([SubjectFile '_' Pfx '_' num2str(idx_rep) '.mat']);
%     system(['cd /home/aglatz/tmp/mineral; ./ps2pdfAll.sh; ' ...
%             'find . -name T2swHypo_mask.nii.gz -execdir pwd \; ' ...
%             '| ./copy_to_host.sh -s brain@dcn060062.dcn.ed.ac.uk:22 ' ...
%             'host/Desktop/Segmentation/AD25_CV_' num2str(idx_rep) '/ ' ...
%             'T2swHypo_mask.nii.gz T2swHypoT1wHypo_mask.nii.gz T2swHypoT1wHyper_mask.nii.gz '...
%             'NormTis_mask.nii.gz NonTis_mask.nii.gz class.pdf Ret.mat']);
end