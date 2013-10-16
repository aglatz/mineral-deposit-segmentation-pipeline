function [H] = print_appvsastructure(Ret, type, fvalid)
N_subject = length(Ret);
N_roi = length(Ret(1).Lab_roi);
% Appearance per structure
Mat_A = NaN(N_subject, N_roi);
for idx_roi = 1:N_roi
    Tmp1 = [Ret.(type)];
    Tmp1 = Tmp1(idx_roi:N_roi:length(Tmp1));
    A = [Tmp1.A_feat];
    M = [Tmp1.I_feat_fvalid] > fvalid;
    Mat_A(M, idx_roi) = A(M);
end
% Calculate occurance per structure in percent
Mat = zeros(3, N_roi);
for idx = 1:N_roi
    Mat(2, idx) = sum(sum(Mat_A(:, idx) == 0));
    Mat(3, idx) = sum(sum(Mat_A(:, idx) == 2));
    Mat(1, idx) = sum(sum(Mat_A(:, idx) == 1));
end
Mat = Mat ./ repmat(sum(Mat, 1), 3, 1) .* 100;
fprintf(['- ' type ' ----------------------------------------------------\n']);
fprintf('Mat_A\n');
(N_subject*ones(1, N_roi) - sum(isnan(Mat_A), 1)) ./ (N_subject*ones(1, N_roi))
fprintf(['---------------------------------------------------------------\n']);
fprintf('Mat - hypo iso hyper\n');
Mat
fprintf(['===============================================================\n']);
H = figure;
HB = bar(1:N_roi, Mat', 'stacked');
axis([0 N_roi+1 0 130]);
set(gca, 'XTick', 1:N_roi);
set(gca, 'XtickLabel', {'Caudate', 'Putamen', 'Globus Pallidus', 'Internal Capsule'});
set(gca, 'YTick', 0:20:100);
set(gca, 'YTickLabel', {'0%', '20%', '40%', '60%', '80%', '100%'});
%set(gca, 'FontSize', 12);
%xlabel('\bf Anatomical structure');
ylabel(['\bf Appearance on ' type ]);
legend_str = cell(1, 3);
legend_str{1, 1} = '\bf hypointense';
legend_str{1, 2} = '\bf isointense';
legend_str{1, 3} = '\bf hyperintense';
%gridLegend(HB, 3, legend_str, 'location', 'southoutside');
legend(legend_str, 'Location', 'NorthWest');
colormap(gray(3));