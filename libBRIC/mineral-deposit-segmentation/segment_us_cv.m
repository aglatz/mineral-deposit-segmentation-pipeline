close all; clear all

SubjectFile = 'subjects_98';
[ndata, text, raw] = xlsread([SubjectFile '.xls']); clear ndata text;
Subjects = raw(:, size(raw, 2));
N_total = length(Subjects);

for idx_rep = 1:3
    cvp = cvpartition(N_total, 'k', 10);
    out = struct;
    RoiLabelTable = {[13, 11, 12, 14]};
    for idx_testset = 1:cvp.NumTestSets
        % Train
        idx_train = cvp.training(idx_testset);
        Subjects_train = Subjects(idx_train);
        SubjectFile_traintest = [SubjectFile '_' num2str(idx_testset)];
        save_xls(SubjectFile_traintest, Subjects_train);
        [~, ~, Over] = segment_us_mp([SubjectFile_traintest '.xls'], ...
                                        RoiLabelTable, 6, 'ThreshFactor', [1 0]);
        M = sum(isnan(Over), 2)==0;
        Ret = mcdregres(Over(M, 1), Over(M,2), 'plots', 0);
        out(idx_testset).P = [Ret.slope Ret.int];

        % Test
        idx_test = cvp.test(idx_testset);
        Subjects_test = Subjects(idx_test);
        save_xls(SubjectFile_traintest, Subjects_test);
        [out(idx_testset).Ret, out(idx_testset).Subjects] = ...
                        segment_us_mp([SubjectFile_traintest '.xls'], ...
                        RoiLabelTable, 6, 'ThreshFactor', out(idx_testset).P);
    end
    save([SubjectFile '_ad25_cv_' num2str(idx_rep) '.mat']);
    system(['cd /home/aglatz/tmp/mineral; ./ps2pdfAll.sh; ' ...
            'find . -name T2swHypo_mask.nii.gz -execdir pwd \; ' ...
            '| ./copy_to_host.sh -s brain@dcn060062.dcn.ed.ac.uk:22 ' ...
            'host/Desktop/Segmentation/AD25_CV_' num2str(idx_rep) '/ ' ...
            'T2swHypo_mask.nii.gz T2swHypoT1wHypo_mask.nii.gz T2swHypoT1wHyper_mask.nii.gz '...
            'NormTis_mask.nii.gz NonTis_mask.nii.gz class.pdf Ret.mat']);
end