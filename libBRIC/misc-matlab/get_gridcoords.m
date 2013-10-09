function [X_old, Y_old, Z_old, X_new, Y_new, Z_new] = get_gridcoords(S_old, F_old, F_new)
% Returns the old and the new lattice coordinates required by interp3().
% INPUTS: S_old - the old volume
%         F_old - the old voxel size
%         F_new - the new voxel size
% OUTPUTS: X_old, Y_old, Z_old - the old lattice coordinates
%          X_new, Y_new, Z_new - the new lattice coordinates
%
% Old lattice coordinates
X_old = 1:F_old(1):(size(S_old, 1)*F_old(1));
Y_old = 1:F_old(2):(size(S_old, 2)*F_old(2));
Z_old = 1:F_old(3):(size(S_old, 3)*F_old(3));
% New lattice coordinates
X_new = 1:F_new(1):(size(S_old, 1)*F_old(1));
Y_new = 1:F_new(2):(size(S_old, 2)*F_old(2));
Z_new = 1:F_new(3):(size(S_old, 3)*F_old(3));
