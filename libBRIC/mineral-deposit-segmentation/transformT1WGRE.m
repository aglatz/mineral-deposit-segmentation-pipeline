% This script processes phantom data. It can remove the appearance of CaAlg
% balls from T1w and replace them with intensities from surrounding agar.
% Then balls are segmented on the co-registered T1w and T2*w volumes and
% some stats are calculated.
% TODO: Needs more documentation.

close all; clear all;

Subject = '/home/aglatz/tmp/mineral-deposit-segmentation-pipeline/libBRIC/qMRI/test/21297/T1W';
% Subject = '/home/aglatz/tmp/mineral-deposit-segmentation-pipeline/libBRIC/qMRI/test/21293/R2s/52716';

% S_t1w = load_series([Subject '/T1W_brain_restore_orig'], []);
% S_gre = load_series([Subject '/GRE_brain_restore_orig'], []);
% 
% save_series([Subject '/T1W_brain_restore_orig'], [Subject '/T1W_brain_restore'], S_t1w, []);
% save_series([Subject '/GRE_brain_restore_orig'], [Subject '/GRE_brain_restore'], S_gre, []);


RoiLabelTable = {1, 2, 4, 8, 16, 32, 64};
% % RoiLabelTable = {1, 2, 3, 4, 5, 6, 7};
% segment_us_single(Subject, RoiLabelTable, 'class', 1, [1 0], true);

% S_roi = load_series([Subject '/RO_mask'], []);
% S_norm = load_series([Subject '/NormTis_mask'], []);
% 
% Lab = cell2mat(RoiLabelTable);
% S_t1w_new = S_t1w;
% for idx_lab = 1:length(Lab)
%     SM_roi = S_roi == Lab(idx_lab);
%     N_roi = sum(SM_roi(:));
%     SM_norm = S_norm == Lab(idx_lab);
%     S_t1w_new(SM_roi) = randsample(S_t1w(SM_norm), N_roi, true);
% end

% S_t1w_new = 400 + 50*randn(size(S_t1w));
% save_series([Subject '/T1W_brain_restore_orig'], [Subject '/T1W_brain_restore'], S_t1w_new, []);

[~, CC] = segment_us_single(Subject, RoiLabelTable, 'class', 1, [1 0], true, 0.3);

S_roi = load_series([Subject '/RO_mask.nii.gz'], []);
SM_hypo = logical(load_series([Subject '/T2swHypo_mask.nii.gz'], []));
NII = load_series([Subject '/RO_mask.nii.gz'], 0);
F = NII.hdr.dime.pixdim(2:4);
Dim = 3;
L_hypo = conncomp_init(SM_hypo, Dim);
SM_hypo_ref = logical(load_series([Subject '/T2swHypo_ref_mask.nii.gz'], []));
L_hypo_ref = conncomp_init(SM_hypo_ref, Dim);

S_t1w = load_series([Subject '/T1W_brain_restore_orig'], []);

Lab = cell2mat(RoiLabelTable);
out = NaN(length(Lab), 8);
VolAll = cell(length(Lab), 1);
It1wAll = cell(length(Lab), 1);
for idx_lab = 1:length(Lab)
    SM_roi = S_roi == Lab(idx_lab);
    Lab_hypo = unique(L_hypo(SM_roi & SM_hypo_ref)); Lab_hypo = Lab_hypo(Lab_hypo > 0);
    Lab_hypo_ref = unique(L_hypo_ref(SM_roi)); Lab_hypo_ref = Lab_hypo_ref(Lab_hypo_ref > 0);
    Lab_hypo_total = unique(L_hypo(SM_roi)); Lab_hypo_total = Lab_hypo_total(Lab_hypo_total > 0);
    Vol = NaN(1, length(Lab_hypo));
    It1w = NaN(1, length(Lab_hypo));
    for idx_lab_cc = 1:length(Lab_hypo)
        SM_tmp = L_hypo == Lab_hypo(idx_lab_cc);
        Vol(idx_lab_cc) = get_volume(SM_tmp, F);
        It1w(idx_lab_cc) = median(S_t1w(SM_tmp));
    end
    if isempty(Vol)
        out(idx_lab, 1) = Lab(idx_lab);
        VolAll{idx_lab} = NaN;
        It1wAll{idx_lab} = NaN;
    else
        out(idx_lab, :) = [Lab(idx_lab) ...
                           length(Lab_hypo_total) length(Lab_hypo) length(Lab_hypo_ref) ...
                           median(Vol) iqr(Vol) mean(Vol) std(Vol)];
        VolAll{idx_lab} = Vol;
        It1wAll{idx_lab} = It1w;
    end
end

save([Subject '/out.mat'], 'CC', 'VolAll', 'out');
csvwrite([Subject '/out.csv'], out);
fprintf('Label Total TP Ref median iqr mean std\n');
out

% Old
%
% SM_roi=logical(load_series([Subject '/RO_mask_AE'], []));
% 
% [I_pc, I_mean, I_sd, V, D, C, C_orig] = pcomp_find([GRE(SM_roi) T1W(SM_roi)]);
% 
% figure;
% scatter(GRE(SM_roi), T1W(SM_roi));
% hold on;
% pcomp_plot(I_mean, I_sd, 'r');
% 
% figure;
% scatter(I_pc(:,1), I_pc(:, 2));
% 
% I = double([GRE(:) T1W(:)]);
% I_tft = ((I - repmat(I_mean, size(I, 1), 1)) * V) / sqrt(D);
% 
% max(I_tft)
% min(I_tft)
% 
% I_tft1=(I_tft+19)/40*4000;
% 
% max(I_tft1)
% min(I_tft1)
% 
% GRE(1:end) = int16(I_tft1(:,1));
% T1W(1:end) = int16(I_tft1(:,2));
