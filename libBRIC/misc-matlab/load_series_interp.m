function S = load_series_interp(fname, slices, method, F)
% Loads 3D volume and interpolates it to the given voxel size.
% Inputs: fname - input volume
%         slices - slices that should be loaded from input volume
%         method - either 'fft' or interpolation method of 'interp3'
%         F - 3 element vector with voxel size
%
S_orig_hdr = load_series(fname, 0);
S_orig_data = load_series(fname, []);
S_orig_type = class(S_orig_data);

dx = S_orig_hdr.hdr.dime.pixdim(2);
dy = S_orig_hdr.hdr.dime.pixdim(3);
dz = S_orig_hdr.hdr.dime.pixdim(4);
fx_inv = get_f_inv(F(1), dx);
fy_inv = get_f_inv(F(2), dy);
fz_inv = get_f_inv(F(3), dz);
X_old = 1:dx*fx_inv:ceil(S_orig_hdr.hdr.dime.dim(2)*dx*fx_inv);
Y_old = 1:dy*fy_inv:ceil(S_orig_hdr.hdr.dime.dim(3)*dy*fy_inv);
Z_old = 1:dz*fz_inv:ceil(S_orig_hdr.hdr.dime.dim(4)*dz*fz_inv);
Z_old_range = [min(Z_old(slices)) max(Z_old(slices))];
X_new = 1:ceil(S_orig_hdr.hdr.dime.dim(2)*dx*fx_inv);
Y_new = 1:ceil(S_orig_hdr.hdr.dime.dim(3)*dy*fy_inv);
Z_new = 1:ceil(S_orig_hdr.hdr.dime.dim(4)*dz*fz_inv);
M_z_new_range = Z_new >= Z_old_range(1) & Z_new <= Z_old_range(2);

switch method
    case 'fft'
        S_tmp = fftInterpolate( double(S_orig_data), ...
                                [length(X_new) length(Y_new) length(Z_new)]);
        S = cast(S_tmp(:, :, M_z_new_range), S_orig_type);
    otherwise
        if isempty(method)
            method = 'linear';
        end
        [Xo, Yo, Zo] = meshgrid(X_old, Y_old, Z_old(slices));
        [Xn, Yn, Zn] = meshgrid(X_new, Y_new, Z_new(M_z_new_range));
        S = cast(interp3(Xo, Yo, Zo, ...
                         double(S_orig_data(:, :, slices)), ...
                         Xn, Yn, Zn, method), S_orig_type);
end


%-------------------------------------------------------------------------
function f_inv = get_f_inv(F, d)
if F == 0
    f_inv = 1/d;
else
    f_inv = 1/F;
end
             
