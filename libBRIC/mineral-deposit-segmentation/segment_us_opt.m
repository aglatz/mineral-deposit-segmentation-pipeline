% This script first estimates the thresholds that lead to the best
% (or at least better) Jaccard indices, which are then used to
% refine the subject thresholds in a 10-fold cross-validation.
%

close all; clear all

SubjectFile = 'subjects_98';
Pfx = 'ad_t1w_mineral_stdint_rep1';
AdaptiveFlag = true;
RoiLabelTable = {[13, 11, 12, 14]};
N_cpus = 6;

% CV
[~, ~, raw] = xlsread([SubjectFile '.xls']);
SubjectsAll = raw(:, size(raw, 2));
N_rep = 1;
N_total = length(SubjectsAll);
for idx_rep = 1:N_rep
    cvp = cvpartition(N_total, 'k', 10);
    out = struct;
    for idx_testset = 1:cvp.NumTestSets
        % Train
        idx_train = cvp.training(idx_testset);
        Subjects_train = SubjectsAll(idx_train);
        SubjectFile_traintest = [SubjectFile '_' Pfx '_' num2str(idx_testset)];
        save_xls(SubjectFile_traintest, {'%s'}, Subjects_train);
       
        DirName = num2str(idx_testset);
        mkdir(DirName);
        IntvarP = 0:0.1:1.5;
        J_mean = zeros(length(IntvarP), 1);
        out(idx_testset).J = cell(length(IntvarP), 1);
        for idx_intvarp = 1:length(IntvarP)
            [Ret, Subjects] = segment_us_mp(  [SubjectFile_traintest '.xls'], ...
                                              RoiLabelTable, N_cpus, ...
                                              'ThreshFactor', [1 0], ...
                                              'AdaptiveFlag', AdaptiveFlag, ...
                                              'SaveMaskFlag', false, ...
                                              'N_gre', 1, ...
                                              'CNR_thr', 0, ...
                                              'phypo_thr', 0.1, ...
                                              'intvar_thr', IntvarP(idx_intvarp));
            MatName = [DirName '/hu_' num2str(IntvarP(idx_intvarp)) '_' Pfx '.mat'];
            save(MatName, 'Ret', 'Subjects');
            [~, J] = load_matdata(MatName, 0);
            out(idx_testset).J{idx_intvarp} = J;
            J_mean(idx_intvarp) = quantile(J, .5);
        end
        
        % Select best
        idx = find(J_mean == max(J_mean));
        [idx, len] = find_largemax(idx, '');
        if mod(len, 2) == 0 % see example of find_largemax() for behaviour
                            % when the number of indices is even
            idx = [idx-1, idx];
        end
        out(idx_testset).IntvarPOpt = mean(IntvarP(idx));
        
        % Test
        idx_test = cvp.test(idx_testset);
        Subjects_test = SubjectsAll(idx_test);
        SubjectFile_traintest = [SubjectFile '_' Pfx '_' num2str(idx_testset)];
        save_xls(SubjectFile_traintest, {'%s'}, Subjects_test);
        
        [out(idx_testset).Ret, out(idx_testset).Subjects] = ...
                        segment_us_mp([SubjectFile_traintest '.xls'], ...
                                      RoiLabelTable, N_cpus, ...
                                      'ThreshFactor', [1 0], ...
                                      'ReportName', 'class', ...
                                      'AdaptiveFlag', AdaptiveFlag, ...
                                      'SaveMaskFlag', true, ...
                                      'N_gre', 1, ...
                                      'CNR_thr', 0, ...
                                      'phypo_thr', 0.1, ...
                                      'intvar_thr', out(idx_testset).IntvarPOpt);
    end
end

MatFileName = [SubjectFile '_' Pfx '.mat'];
save(MatFileName, 'out');