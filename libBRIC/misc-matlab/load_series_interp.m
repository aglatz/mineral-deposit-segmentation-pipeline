function S = load_series_interp(fname, slices, method, F)
% Loads 3D volume and interpolates it to the given voxel size.
% Inputs: fname - input volume
%         slices - slices that should be loaded from input volume
%         method - interpolation method (see help interp3)
%         F - 3 element vector with voxel size
%
S_orig_hdr = load_series(fname, 0);
S_orig_data = load_series(fname, slices);

if isempty(slices)
    z_max = S_orig_hdr.hdr.dime.dim(4);
else
    z_max = length(slices);
end
dx = S_orig_hdr.hdr.dime.pixdim(2);
dy = S_orig_hdr.hdr.dime.pixdim(3);
dz = S_orig_hdr.hdr.dime.pixdim(4);
fx_inv = get_f_inv(F(1), dx);
fy_inv = get_f_inv(F(2), dy);
fz_inv = get_f_inv(F(3), dz);
[Xo, Yo, Zo] = meshgrid(   1:dx*fx_inv:ceil(S_orig_hdr.hdr.dime.dim(2)*dx*fx_inv), ...
                           1:dy*fy_inv:ceil(S_orig_hdr.hdr.dime.dim(3)*dy*fy_inv), ...
                           1:dz*fz_inv:ceil(z_max*dz*fz_inv));
[Xn, Yn, Zn] = meshgrid(   1:ceil(S_orig_hdr.hdr.dime.dim(2)*dx*fx_inv), ...
                           1:ceil(S_orig_hdr.hdr.dime.dim(3)*dy*fy_inv), ...
                           1:ceil(z_max*dz*fz_inv));

if isempty(method)
    method = 'linear';
end
type = class(S_orig_data);
S = cast(interp3(Xo, Yo, Zo, double(S_orig_data), Xn, Yn, Zn, method), type);


%-------------------------------------------------------------------------
function f_inv = get_f_inv(F, d)
if F == 0
    f_inv = 1/d;
else
    f_inv = 1/F;
end
             
