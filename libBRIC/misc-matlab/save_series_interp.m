function [] = save_series_interp(fname, new_name, S_orig_data, slices, method, F)
% Interpolates the given volume to the given voxel size and saves it.
% Inputs: fname - template volume
%         new_name - new volume name
%         S_orig_data - Data to save
%         slices - slices that were loaded from original volume or empty
%         method - interpolation method of 'interp3'
%         F - 3 element vector with voxel size
%
S_orig_type = class(S_orig_data);

[X_old, Y_old, Z_old, X_new, Y_new, Z_new] = get_gridcoords(fname, F);
Z_old_range = [min(Z_old(slices)) max(Z_old(slices))];
M_z_new_range = Z_new >= Z_old_range(1) & Z_new <= Z_old_range(2);

[Xo, Yo, Zo] = meshgrid(X_old, Y_old, Z_old(slices));
[Xn, Yn, Zn] = meshgrid(X_new, Y_new, Z_new(M_z_new_range));
if isempty(method)
    method = 'linear';
end
S = cast(interp3(Xn, Yn, Zn, double(S_orig_data), Xo, Yo, Zo, method), S_orig_type);
save_series(fname, new_name, S, slices);
