function [H] = print_appvsastructure(Ret, type)
N_subject = length(Ret);
N_roi = 4; %length(Ret(1).Lab_roi);
Tmp = [Ret.(type)];

% Calculate occurance per structure in percent
Tmp = reshape([Tmp.T_roi_fe], length(Ret(1).Lab_roi), N_subject)';
Mat = zeros(3, N_roi);
for idx = 1:N_roi
    Mat(2, idx) = sum(sum(Tmp(:, [idx idx+4]) == 0));
    Mat(3, idx) = sum(sum(Tmp(:, [idx idx+4]) == 2));
    Mat(1, idx) = sum(sum(Tmp(:, [idx idx+4]) == 1));
end

Mat = Mat ./ repmat(sum(Mat, 1), 3, 1) .* 100;
fprintf(['- ' type ' ----------------------------------------------------\n']);
fprintf('Mat - hypo iso hyper\n');
Mat
fprintf(['===============================================================\n']);
H = figure;
HB = bar(1:N_roi, Mat', 'stacked');
axis([0 N_roi+1 0 130]);
set(gca, 'XTick', 1:N_roi);
set(gca, 'XtickLabel', {'Caudate', 'Putamen', 'Globus Pallidus', 'Internal Capsule'}); %{'CL', 'PL', 'GL', 'IL', 'CR', 'PR', 'GR', 'IR'});
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