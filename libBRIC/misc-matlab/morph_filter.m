function [SM_oli] = morph_filter(fname, slices, SM_oli)
% Morphological filtering of mask
%
fprintf('Before: %d ...', sum(SM_oli(:)));
L = conncomp_init(SM_oli, 3);
S_gre_hdr = load_series(fname, 0);
dx = S_gre_hdr.hdr.dime.pixdim(3);
dy = S_gre_hdr.hdr.dime.pixdim(2);
dz = S_gre_hdr.hdr.dime.pixdim(4);
S_gre = double(load_series(fname, slices));
CC = conncomp_list(L, L, SM_oli, S_gre, [dx dy dz]);
if CC(1).lab > -1
    M = [CC.cont] < -0.001; % & ...
        %(([CC.nsl] == 1 & [CC.vol] >= 2) | ([CC.nsl] > 1 & [CC.vol] >= 4));
    Lab = [CC(M).lab]';
    N_lab = length(Lab);
    SM_oli = false(size(SM_oli));
    for idx_lab = 1:N_lab
        SM_oli = SM_oli | L == Lab(idx_lab);
    end
end
fprintf('After: %d.\n', sum(SM_oli(:)));