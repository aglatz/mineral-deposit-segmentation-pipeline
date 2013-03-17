function [] = pcomp_plot(I_mean, I_sd, col)
% Plots the principal components.
% INPUTS: I_mean - the mean of the bivariate distribution
%         I_sd - the standard deviation calculated from
%                the covariance matrix in the coordinate
%                system given by the principal components
%         col - color
%
hold on;
for i = 1:2
    Comp = [I_mean(:) I_mean(:) + I_sd(:, i)];
    plot(Comp(1, :), Comp(2, :), 'Color', col);
end
