function [H] = print_intvsastructure(Ret, type, Q_idx, fvalid, varargin)
N_subject = length(Ret);
N_roi = length(Ret(1).Lab_roi);
Mat_ntis = NaN(N_subject, N_roi);
Mat_feat = NaN(N_subject, N_roi);
for idx_roi = 1:N_roi
    % Normal appearing tissue intensities
    Tmp1 = [Ret.(type)];
    Tmp1 = Tmp1(idx_roi:N_roi:length(Tmp1));
    Tmp2 = reshape([Tmp1.I_ntis], length(Tmp1(1).I_ntis), length(Tmp1))';
    M = [Tmp1.I_ntis_fvalid] > fvalid;
    Mat_ntis(M, idx_roi) = Tmp2(M, Q_idx);
    % Feature tissue intensities
    Tmp2 = reshape([Tmp1.I_feat], length(Tmp1(1).I_feat), length(Tmp1))';
    M = [Tmp1.I_feat_fvalid] > fvalid;
    Mat_feat(M, idx_roi) = Tmp2(M, Q_idx);
end
if ~isempty(varargin)
    Mat_ref = repmat(Mat_ntis(:, varargin{1}), 1, size(Mat_ntis, 2));
    Mat_ntis = Mat_ntis ./ Mat_ref;
    Mat_feat = Mat_feat ./ Mat_ref;
end
% Average values
fprintf(['- ' type ' ----------------------------------------------------\n']);
fprintf('Mat_ntis\n');
(N_subject*ones(1, N_roi) - sum(isnan(Mat_ntis), 1)) ./ (N_subject*ones(1, N_roi))
quantile(Mat_ntis, [.25 .5 .75])
fprintf(['---------------------------------------------------------------\n']);
fprintf('Mat_feat\n');
(N_subject*ones(1, N_roi) - sum(isnan(Mat_feat), 1)) ./ (N_subject*ones(1, N_roi))
quantile(Mat_feat, [.25 .5 .75])
fprintf(['---------------------------------------------------------------\n']);
% Correlation
[R, P] = corr(quantile(Mat_ntis, [.5])', quantile(Mat_feat, [.5])', 'type', 'Spearman')
fprintf(['---------------------------------------------------------------\n']);
for idx=1:N_roi
    fprintf('%d:', idx);
    Col_idx = idx;
    Tmp_ntis= Mat_ntis(:, Col_idx);
    Tmp_feat = Mat_feat(:, Col_idx);
    Tmp = quantile((Tmp_ntis(:)-Tmp_feat(:))./Tmp_ntis(:), [.25 .5 .75]);
    [Tmp(2) (Tmp(3)-Tmp(1))/2]
end
quantile((Mat_ntis(:)-Mat_feat(:))./Mat_ntis(:), [.25 .5 .75])
fprintf(['---------------------------------------------------------------\n']);
for idx=1:N_roi
    fprintf('%d:', idx);
    Col_idx = idx;
    Tmp_feat = Mat_feat(:, Col_idx);
    quantile(Tmp_feat(:), [.25 .5 .75])
end
quantile(Mat_feat(:), [.25 .5 .75])
fprintf(['---------------------------------------------------------------\n']);
Col_idx_ref = 1; % caudate
Tmp_ref = Mat_ntis(:, Col_idx_ref);
for idx=1:N_roi-1
    fprintf('%d:', idx);
    Col_idx = idx+1;
    Tmp_roi = Mat_ntis(:, Col_idx);
    Tmp = quantile((Tmp_ref(:)-Tmp_roi(:))./Tmp_ref(:), [.25 .5 .75]);
    [Tmp(2) (Tmp(3)-Tmp(1))/2]
end
fprintf(['===============================================================\n']);
Mat = [Mat_ntis; Mat_feat];
G1 = [ones(N_subject, N_roi); 2*ones(N_subject, N_roi)];
G2 = repmat([1 2 3 4], 2*N_subject, 1);
H = figure;
CM = gray(3);
hl = boxplot(Mat(:), {G2(:), G1(:)}, 'colorgroup', G1(:), ...
             'factorgap', 5, 'factorseparator', 1, ...
             'plotstyle', 'compact', 'labelorientation', 'horizontal', ...
             'labels', {'Caudate', '', 'Putamen', '', 'Globus', 'Pallidus', 'Internal', 'Capsule'}, ...
             'labelverbosity', 'all', ...
             'colors', CM(1:2, :));
ylabel(['\bf Standardized ' type ' intensities in arbitrary units']);
Obj = findobj(gca,'Tag','Box');
legend(Obj(1:2), 'Mineral deposits', 'Normal tissue', 'Location', 'Best');
