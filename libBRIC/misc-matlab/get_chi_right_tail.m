function [Area] = get_chi_right_tail(I, I_mean, C)
% Return Area under right tail of chi2 distribution
% INPUTS I - Matrix of T1- and T2*-weighted signal intensities
%        I_mean - Vector with normal-appearing tissue intensities
%        C - Covariance matrix of normal-appearing tissue intensities
% OUTPUT Area - Area under chi2 distribution
%
I_new = I - repmat(I_mean(:), 1, size(I, 2));
Area = ((C(2,2).*I_new(1, :)-C(1,2).*I_new(2, :)).^2./det(C) + I_new(2, :).^2)./C(2,2);
