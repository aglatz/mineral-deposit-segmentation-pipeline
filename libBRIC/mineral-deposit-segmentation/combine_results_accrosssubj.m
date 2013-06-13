close all; clear all

Name = 'subjects_98';

[J, D, V, ~, ~, FPC, FPV, ~] = load_matdata([Name '_fi25.mat']);
N_ex = 8;
idx = 1;
J_all = NaN(length(J), N_ex);
J_all(:, idx) = J;
D_all = NaN(length(D), N_ex);
D_all(:, idx) = D;
V_all = NaN(length(V), N_ex);
V_all(1:length(FPV), idx) = FPV;

MatNames = {[''], ...
            [Name '_fi25_cv_1.mat'], ...
            [Name '_fi25_cv_2.mat'], ...
            [Name '_fi25_cv_3.mat'], ...
            [Name '_ad25.mat'], ...
            [Name '_ad25_cv_1.mat'], ...
            [Name '_ad25_cv_2.mat'], ...
            [Name '_ad25_cv_3.mat']};

for idx = 2:N_ex
    [J, D, V, ~, ~, FPC, FPV, ~] = load_matdata(MatNames{idx});
    J_all(:, idx) = J;
    D_all(:, idx) = D;
    V_all(1:length(FPV), idx) = FPV;
end

G1 = [1 2 2 2 3 4 4 4];
G2 = 1:8;
CM = gray(5);
figure;
subplot(311);
boxplot(J_all, {G1(:) G2(:)}, 'colorgroup', G1(:), ...
             'factorgap', 5, 'factorseparator', 1, ...
             'plotstyle', 'compact', 'labelorientation', ...
             'horizontal',  ...
             'colors', CM([2 3 4 1], :));
ylabel('\bf Jaccard index');
ylim([0 1]);
subplot(312);
boxplot(D_all, {G1(:) G2(:)}, 'colorgroup', G1(:), ...
             'factorgap', 5, 'factorseparator', 1, ...
             'plotstyle', 'compact', 'labelorientation', ...
             'horizontal',  ...
             'colors', CM([2 3 4 1], :));
ylabel('\bf Dice coefficient');
ylim([0 1]);
subplot(313);
boxplot(V_all, {G1(:) G2(:)}, 'colorgroup', G1(:), ...
             'factorgap', 5, 'factorseparator', 1, ...
             'plotstyle', 'compact', 'labelorientation', ...
             'horizontal',  ...
             'colors', CM([2 3 4 1], :));
ylabel('\bf FP volume in mm^3');
set(gca, 'YScale', 'log');
set(gcf, 'color', 'white');

