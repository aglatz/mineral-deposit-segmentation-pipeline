function [out1, out2] = transformT1WGRE(Subject, adapt_flag, te, intstd_thr)

% This script processes phantom data. It can remove the appearance of CaAlg
% balls from T1w and replace them with intensities from surrounding agar.
% Then balls are segmented on the co-registered T1w and T2*w volumes and
% some stats are calculated.
% TODO: Needs more documentation.

% close all; clear all;

% Subject = '/home/aglatz/tmp/mineral-deposit-segmentation-pipeline/libBRIC/qMRI/test/21297/T1W';
% Subject = '/home/aglatz/tmp/mineral-deposit-segmentation-pipeline/libBRIC/qMRI/test/21293/R2s/52716';

% S_t1w = load_series([Subject '/T1W_brain_restore_orig'], []);
% S_gre = load_series([Subject '/GRE_brain_restore_orig'], []);
% 
% save_series([Subject '/T1W_brain_restore_orig'], [Subject '/T1W_brain_restore'], S_t1w, []);
% save_series([Subject '/GRE_brain_restore_orig'], [Subject '/GRE_brain_restore'], S_gre, []);


RoiLabelTable = {1, 2, 3, 11, 12, 13, 21, 22, 23};
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

Ret = segment_us_single(Subject, RoiLabelTable, 'class', 1, [1 0], adapt_flag, te, 0, 0.1, intstd_thr);
I_thr = [Ret.I_thr{:}]; % first row: gre thr

S_roi = load_series([Subject '/RO_mask.nii.gz'], []);
S_norm = load_series([Subject '/NormTis_mask.nii.gz'], []);
SM_hypo = logical(load_series([Subject '/T2swHypo_mask.nii.gz'], []));
NII = load_series([Subject '/RO_mask.nii.gz'], 0);
F = NII.hdr.dime.pixdim(2:4);
Dim = 3;
% L_hypo = conncomp_init(SM_hypo, Dim);
SM_hypo_ref = logical(load_series([Subject '/FE_roi_mask.nii.gz'], []));
L_hypo_ref = conncomp_init(SM_hypo_ref, Dim);
%save_series([Subject '/FE_roi_mask.nii.gz'], [Subject '/FE_roi_mask_cc.nii.gz'], L_hypo_ref, []);

S_gre = single(load_series([Subject '/GRE_brain_restore'], []));
fd = fopen([Subject '/GRE_brain_restore_mask.txt'], 'r');
Noise_gre = single(fscanf(fd, '%d'));
fclose(fd);

S_t1w = single(load_series([Subject '/T1w_brain_restore'], []));
fd = fopen([Subject '/T1w_brain_restore_mask.txt'], 'r');
Noise_t1w = single(fscanf(fd, '%d'));
fclose(fd);

Lab = [1,2,3,11,12,13,21,22,23]; % dont change order
out1 = NaN(length(Lab), 5);
out2 = NaN(length(Lab), 7, 4);
for idx_lab = 1:length(Lab)
	SM_roi = S_roi == Lab(idx_lab);
    SM_norm = S_norm == Lab(idx_lab);
    out1(idx_lab, 1) = Lab(idx_lab);
    Lab_hypo = unique(L_hypo_ref(SM_roi & SM_hypo)); Lab_hypo = Lab_hypo(Lab_hypo > 0);
    out1(idx_lab, 2) = length(Lab_hypo);
    Lab_hypo_ref = unique(L_hypo_ref(SM_roi)); Lab_hypo_ref = Lab_hypo_ref(Lab_hypo_ref > 0);
    out1(idx_lab, 3) = length(Lab_hypo_ref);
    for idx2_lab = 1:7 %length(Lab_hypo_ref)
        SM_tmp = L_hypo_ref == Lab_hypo_ref(idx2_lab);
        % position
        [~, ~, z] = ind2sub(size(SM_tmp), find(SM_tmp));
        out2(idx_lab, idx2_lab, 1) = median(z);
        % gre contrast
        out2(idx_lab, idx2_lab, 2) = abs( median(S_gre(SM_tmp)) - ...
                                          median(S_gre(SM_norm))) / Noise_gre;
        if (~sum(Lab_hypo == Lab_hypo_ref(idx2_lab)))
            out2(idx_lab, idx2_lab, 2) = out2(idx_lab, idx2_lab, 2)*-1;
        end
        % t1w contrast
        out2(idx_lab, idx2_lab, 3) = abs( median(S_t1w(SM_tmp)) - ...
                                          median(S_t1w(SM_norm))) / Noise_t1w;
        if (~sum(Lab_hypo == Lab_hypo_ref(idx2_lab)))
            out2(idx_lab, idx2_lab, 3) = out2(idx_lab, idx2_lab, 3)*-1;
        end
        % size
        out2(idx_lab, idx2_lab, 4) = get_volume(SM_tmp, F);
    end
    out1(idx_lab, 5) = abs(I_thr(1, idx_lab) - median(S_gre(SM_norm))) / Noise_gre;
end

for idx_lab = 1:3:length(Lab)
    SM_roi = S_roi == Lab(idx_lab);
    for idx1 = 1:2
        SM_roi = SM_roi | S_roi == Lab(idx_lab+idx1);
    end
    Ret = validate_raw(SM_roi & SM_hypo, SM_roi & SM_hypo_ref, [], F);
    out1(idx_lab, 4) = Ret.Jaccard;    
end

fprintf('Label Selected Ref\n');
Tmp=out1(1:3, :); [quantile(Tmp, .5); iqr(Tmp)/ 1.349]
Tmp=out1(4:6, :); [quantile(Tmp, .5); iqr(Tmp)/ 1.349]
Tmp=out1(7:9, :); [quantile(Tmp, .5); iqr(Tmp)/ 1.349]
csvwrite([Subject '/out1.csv'], out1);
out2
csvwrite([Subject '/out2.csv'], out2);
fprintf('------------------\n');

figure;
subplot(3,1,1);
boxplot(out1(:, [2]),[1 1 1 2 2 2 3 3 3], 'labels', {'N/1200nm+HA', 'N/250nm', 'N/1200nm'});
ylabel('Count');
subplot(3,1,2);
boxplot(out1(:, [4]),[1 1 1 2 2 2 3 3 3], 'labels', {'N/1200nm+HA', 'N/250nm', 'N/1200nm'});
ylabel('Jaccard index');
subplot(3,1,3);
boxplot(out1(:, [5]),[1 1 1 2 2 2 3 3 3], 'labels', {'N/1200nm+HA', 'N/250nm', 'N/1200nm'});
ylabel('Threshold (CNR)');

conc = [7 5 4 3 2 1 0];
for idx1 = 1:3
    figure
    hold on;
    C = hsv(3);
    L = cell(3, 1);
    for idx = 1:3
        Tmp = abs(out2(1:3+(idx-1)*3, :, 2+idx1-1));
        if idx1 == 3
            Tmp = (Tmp*3/4/pi).^(1/3);
        end
        Y = median(Tmp);
        Y_std = iqr(Tmp) ./ 1.349;
        errorbar(conc, Y(1:length(conc)), Y_std(1:length(conc)), 'color', C(idx, :));
        L{idx} = sprintf('%d,', out1(1+(idx-1)*3, 1));
    end
    legend('N/1200nm+HA', 'N/250nm', 'N/1200nm');
%     if idx1 == 1
%         for idx = 1:3
%             Thr = median(out3(1:3+(idx-1)*3));
%             plot(conc([1 end]), [Thr, Thr], 'color', C(idx, :));
%             fprintf('%s: %0.2f\n', L{idx}, Thr);
%         end
%     end
    xlabel('\bf Concentration in mmol/l');
    set(gcf, 'color', 'white');
    if idx1 == 3
        ylabel('\bf Volume in mm^3');
    else
        ylabel('\bf Contrast to noise ratio');
        axis([0     8     0     10]);
    end
end
fprintf('------------------\n');



% S_t1w = load_series([Subject '/T1W_brain_restore'], []);
% 
% Lab = cell2mat(RoiLabelTable);
% out = NaN(length(Lab), 8);
% VolAll = cell(length(Lab), 1);
% It1wAll = cell(length(Lab), 1);
% for idx_lab = 1:length(Lab)
%     SM_roi = S_roi == Lab(idx_lab);
%     Lab_hypo = unique(L_hypo(SM_roi & SM_hypo_ref)); Lab_hypo = Lab_hypo(Lab_hypo > 0);
%     Lab_hypo_ref = unique(L_hypo_ref(SM_roi)); Lab_hypo_ref = Lab_hypo_ref(Lab_hypo_ref > 0);
%     Lab_hypo_total = unique(L_hypo(SM_roi)); Lab_hypo_total = Lab_hypo_total(Lab_hypo_total > 0);
%     Vol = NaN(1, length(Lab_hypo));
%     It1w = NaN(1, length(Lab_hypo));
%     for idx_lab_cc = 1:length(Lab_hypo)
%         SM_tmp = L_hypo == Lab_hypo(idx_lab_cc);
%         Vol(idx_lab_cc) = get_volume(SM_tmp, F);
%         It1w(idx_lab_cc) = median(S_t1w(SM_tmp));
%     end
%     if isempty(Vol)
%         out(idx_lab, 1) = Lab(idx_lab);
%         VolAll{idx_lab} = NaN;
%         It1wAll{idx_lab} = NaN;
%     else
%         out(idx_lab, :) = [Lab(idx_lab) ...
%                            length(Lab_hypo_total) length(Lab_hypo) length(Lab_hypo_ref) ...
%                            median(Vol) iqr(Vol) mean(Vol) std(Vol)];
%         VolAll{idx_lab} = Vol;
%         It1wAll{idx_lab} = It1w;
%     end
% end
% 
% save([Subject '/out.mat'], 'VolAll', 'out');
% csvwrite([Subject '/out.csv'], out);
% fprintf('Label Total TP Ref median iqr mean std\n');
% out

%NN = zeros(5, 1);


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
