function [X_old, Y_old, Z_old, X_new, Y_new, Z_new] = get_gridcoords(fname, F)
% Returns the old and the new lattice coordinates required by interp3().
% INPUTS: fname - the filename of the volume to interpolate
%         F - the interpolation factor (e.g. 2 increase lattice by 200%)
% OUTPUTS: X_old, Y_old, Z_old - the old lattice coordinates
%          X_new, Y_new, Z_new - the new lattice coordinates
%
S_orig_hdr = load_series(fname, 0);
dx = S_orig_hdr.hdr.dime.pixdim(3);         
dy = S_orig_hdr.hdr.dime.pixdim(2);
dz = S_orig_hdr.hdr.dime.pixdim(4);
% Old lattice coordinates
X_old = 1:dx:(S_orig_hdr.hdr.dime.dim(3)*dx);
Y_old = 1:dy:(S_orig_hdr.hdr.dime.dim(2)*dy);
Z_old = 1:dz:(S_orig_hdr.hdr.dime.dim(4)*dz);
% New lattice coordinates
X_new = 1:(dx/F):(S_orig_hdr.hdr.dime.dim(3)*dx);
Y_new = 1:(dy/F):(S_orig_hdr.hdr.dime.dim(2)*dy);
Z_new = 1:(dz/F):(S_orig_hdr.hdr.dime.dim(4)*dz);
