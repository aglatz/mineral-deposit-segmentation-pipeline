function S = load_series_interp(fname, slices, method, F)
% Loads 3D volume and interpolates it to the given voxel size.
% Inputs: fname - input volume
%         slices - slices that should be loaded from input volume
%         method - either 'fft' or interpolation method of 'interp3'
%         F - 3 element vector with voxel size
%
S_orig_data = load_series(fname, []);
S_orig_type = class(S_orig_data);

[X_old, Y_old, Z_old, X_new, Y_new, Z_new] = get_gridcoords(fname, F);
Z_old_range = [min(Z_old(slices)) max(Z_old(slices))];
M_z_new_range = Z_new >= Z_old_range(1) & Z_new <= Z_old_range(2);

switch method
    case 'fft'
        S_tmp = fftInterpolate( double(S_orig_data), ...
                                [length(Y_new) length(X_new) length(Z_new)]);
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
       
