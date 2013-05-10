function [X_old, Y_old, Z_old, X_new, Y_new, Z_new] = get_gridcoords(fname, F)
S_orig_hdr = load_series(fname, 0);
dx = S_orig_hdr.hdr.dime.pixdim(3);         
dy = S_orig_hdr.hdr.dime.pixdim(2);
dz = S_orig_hdr.hdr.dime.pixdim(4);
fx_inv = get_f_inv(F(1), dx);
fy_inv = get_f_inv(F(2), dy);
fz_inv = get_f_inv(F(3), dz);
X_old = 1:dx*fx_inv:ceil(S_orig_hdr.hdr.dime.dim(3)*dx*fx_inv);
Y_old = 1:dy*fy_inv:ceil(S_orig_hdr.hdr.dime.dim(2)*dy*fy_inv);
Z_old = 1:dz*fz_inv:ceil(S_orig_hdr.hdr.dime.dim(4)*dz*fz_inv);
X_new = 1:ceil(S_orig_hdr.hdr.dime.dim(3)*dx*fx_inv);
Y_new = 1:ceil(S_orig_hdr.hdr.dime.dim(2)*dy*fy_inv);
Z_new = 1:ceil(S_orig_hdr.hdr.dime.dim(4)*dz*fz_inv);

%-------------------------------------------------------------------------
function f_inv = get_f_inv(F, d)
if F == 0
    f_inv = 1/d;
else
    f_inv = 1/F;
end