clear all; close all

Mu1 = [0 0]';
Sig1 = [[1 -0.4]; [-0.4 .8]];%diag([1 1]);
N = 1000;
X1 = mvnrnd(Mu1, Sig1, round(85/100*N));

Mu2 = [-4 0]';
Sig2 = diag([1/10 1/10]);
X2 = mvnrnd(Mu2, Sig2, round(15/100*N));

figure;
subplot(211);
scatter(X1(:, 1), X1(:, 2), 10, 'b');
hold on;
scatter(X2(:, 1), X2(:, 2), 10, 'r');

X = [X1; X2];
X = X(randperm(size(X, 1)), :);
[~, I_mean, ~, ~, ~, C, ~] = pcomp_find(X);
MDsqrt = mahalanobis(X, I_mean, 'cov', C);
subplot(212);
[Chi2cest, BinCentre] = ecdf(MDsqrt);
plot(BinCentre, Chi2cest, 'b');
hold on
Chi2c = chi2cdf(BinCentre, 2);
plot(BinCentre, Chi2c, 'r');
cutoff = chi2inv(0.975, 2); % cutoff - fixed sq. MD from chi2 dist
vline(cutoff, 'b', '');
cutoff_final = cutoff; % cutoff_final - adaptive sq. MD cutoff
d = Chi2c - Chi2cest;
dsup = max(d(BinCentre >= cutoff & d >= 0));
if ~isempty(dsup)
    cutoff_adj = BinCentre(find(Chi2cest > (1-dsup), 1));
    vline(cutoff_adj, 'g', '');
    pcrit = (0.24 - 0.003 *2)/sqrt(size(X, 1));
    if dsup > pcrit
        fprintf('Adjust!\n');
        cutoff_final = cutoff_adj;
    end
end

figure;
subplot(311);
scatter(X1(:, 1), X1(:, 2), 10, 'b');
hold on;
scatter(X2(:, 1), X2(:, 2), 10, 'r');
ylabel(sprintf('\\bf Standardized\nT1w signal intensities'));
xlabel('\bf Standardized T2*w signal intensities');

subplot(312);
scatter(X(:, 1), X(:, 2), 10, 'b');
hold on
M = MDsqrt > cutoff;
scatter(X(M, 1), X(M, 2), 10, 'r');
x1_minmax = sqrt(cutoff*C(1,1));
x1 = (-x1_minmax+eps:0.01:x1_minmax+0.01);
x21 = (C(1,2)/C(1,1)).*x1 + sqrt( ((C(1,2)/C(1,1))^2 - C(2,2)/C(1,1)).*(x1.^2) + det(C)*cutoff/C(1,1) );
x22 = (C(1,2)/C(1,1)).*x1 - sqrt( ((C(1,2)/C(1,1))^2 - C(2,2)/C(1,1)).*(x1.^2) + det(C)*cutoff/C(1,1) );
plot(x1 + I_mean(1), x21 + I_mean(2), 'g');
plot(x1 + I_mean(1), x22 + I_mean(2), 'g');
% ellipsplot_mod(I_mean', C, [], 'g', cutoff);
hline(sqrt(cutoff*C(2,2)) + I_mean(2), 'k', '');
hline(-sqrt(cutoff*C(2,2)) + I_mean(2), 'k', '');
vline(-sqrt(cutoff*C(1,1)) + I_mean(1), 'k', num2str(sqrt(cutoff)));
ylabel(sprintf('\\bf Standardized\nT1w signal intensities'));
xlabel('\bf Standardized T2*w signal intensities');

subplot(313);
scatter(X(:, 1), X(:, 2), 10, 'b');
hold on
M = MDsqrt > cutoff_final;
scatter(X(M, 1), X(M, 2), 10, 'r');
hold on;

% PLot ellipse
x1_minmax = sqrt(cutoff_final*C(1,1));
x1 = (-x1_minmax+eps:0.01:x1_minmax+0.01);
x21 = (C(1,2)/C(1,1)).*x1 + sqrt( ((C(1,2)/C(1,1))^2 - C(2,2)/C(1,1)).*(x1.^2) + det(C)*cutoff_final/C(1,1) );
x22 = (C(1,2)/C(1,1)).*x1 - sqrt( ((C(1,2)/C(1,1))^2 - C(2,2)/C(1,1)).*(x1.^2) + det(C)*cutoff_final/C(1,1) );
plot(x1 + I_mean(1), x21 + I_mean(2), 'g');
plot(x1 + I_mean(1), x22 + I_mean(2), 'g');
% ellipsplot_mod(I_mean', C, [], 'g', cutoff_final);

% x1_max = sqrt((cutoff_final*C(1,2).^2)/C(2,2));
% vline(x1_max + I_mean(1), 'g', '');
% vline(-x1_max + I_mean(1), 'g', '');

hline(sqrt(cutoff_final*C(2,2)) + I_mean(2), 'k', '');
hline(-sqrt(cutoff_final*C(2,2)) + I_mean(2), 'k', '');
vline(-sqrt(cutoff_final*C(1,1)) + I_mean(1), 'k', num2str(sqrt(cutoff_final)));
ylabel(sprintf('\\bf Standardized\nT1w signal intensities'));
xlabel('\bf Standardized T2*w signal intensities');






