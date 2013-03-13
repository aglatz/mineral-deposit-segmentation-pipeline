function [I_pc, I_mean, I_sd, V, D, C, C_orig] = pcomp_find(I)
% Finds principal components of bivariate distribution I using
% MCDCOV (LIBRA toolbox; Verboven, S., Chemometr Intell Lab, 2002)
% for estimating the robust mean and covariance matrix.
% INPUTS: I - Nx2 matrix with values of a bivariate distribution
% RETURNS: I_pc - Transformed 'I' using principal component (PC) vectors
%          I_mean - robust mean
%          I_std - robust standard deviations along PCs
%          V - PCs (Eigenvectors of covariance matrix)
%          D - Eigenvalues matrix
%          C - robust covariance matrix after reorientation of PCs so that
%              the point in the direction of the nearest axis of original
%              coordinate system.
%          C_orig - covariance matrix from MCDCOV
%
% To understand what 'reorienting along the nearest axis of original
% coordinate sytem' means run the following code:
% mu = [1 -1];
% f = 0.2; % if f>0.5 the sign of a changes
% figure;
% for a = [-150 -120 -60 -30 30 60 120 150]
% V = [[cosd(a) -sind(a)]; [sind(a) cosd(a)]];
% D = [[(1-f)^2 0]; [0 f^2]];
% Sigma = V*D*V';
% I = mvnrnd(mu, Sigma, 500);
% [I_pc, I_mean, I_sd, V_est, D_est, C] = pcomp_find(I);
% subplot(2, 1, 1);
% scatter(I(:, 1), I(:, 2));
% pcomp_plot(I_mean, I_sd, 'k');
% ellipsplot(I_mean, C, I, 'r', 10.59663);
% axis equal;
% hold off;
% title(sprintf('%0.1f', a));
% xlabel('1st column');
% ylabel('2nd column');
% subplot(2, 1, 2);
% scatter(I_pc(:, 1), I_pc(:, 2));
% pcomp_plot([0 0], [[1 0]; [0 1]], 'k');
% ellipsplot([0 0], [[1 0]; [0 1]], I_pc, 'r', 10.59663);
% axis equal;
% xlabel('1st principal component');
% ylabel('2nd principal component');
% hold off;
% input('');
% end
%
I = double(I);
C_orig = mcdcov(I, 'plots', 0, 'alpha', 0.5);
[V, D] = eig(C_orig.cov);

% We want the principal component first
if D(1, 1) < D(2, 2)
    V(:, [1 2]) = V(:, [2 1]);
    D([1 2],[1 2]) = D([2 1],[2 1]);
end

% Unity vectors of new coordinate system should point in same
% direction as old coordinate system unity vectors.
if sign(V(1, 1)) < 0
    V(:, 1) = -1 .* V(:, 1);
end
if sign(V(2, 2)) < 0
    V(:, 2) = -1 .* V(:, 2);
end

% V must be positive definite
if det(V) < 0
    error('Eigenvector matrix is negative definite!');
end

% Recalculate cov matrix accordingly after rearranging things
C = V*D*V';

% Get mean and standard deviation of components in original space
I_mean = C_orig.center;
I_sd = V*sqrt(D);

% Intensities in standardized PC coordinate system, 
% i.e. mu = [0 0]; Sigma = [[1 0]; [0 1]];
I_pc = ((I - repmat(I_mean, size(I, 1), 1)) * V) / sqrt(D);
