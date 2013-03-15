function [S_lab, Lab] = conncomp_init(SM, dim)
% Calculates the connected components of the given binary mask.
% INPUTS: SM - logical mask
%         dim - specifies the neighbourhood. Possible values:
%               2...2D eight connected neighbourhood
%               3...3D six connected neighbourhood
% RETURNS: S_lab - a matrix the size of 'SM' where each
%                  connected component has its own label
%                  (see: help bwlabeln).
%          Lab - List of unique connected component labels
%
switch dim
    case 2,
        % Just in-plane connectedness because we
        % might have a slice gap or the slice is
        % really thick.
        SM_kern = false(3, 3, 3);
        SM_kern = true(size(SM_kern(:, :, 2)));
    case 3,
        % In 3D it must be 6 so that the compactness
        % and relative anisotropy doesn't get negative.
        SM_kern = 6;
    otherwise,
        error('dim parameter out of range');
end
S_lab = bwlabeln(SM, SM_kern);
Lab_tmp = unique(S_lab)';
Lab = Lab_tmp(2:end); % don't count background

