function [S] = interp_series(S_orig_data, slices, method, F_old, F_new)
% Interpolates 3D volume to the given voxel size.
% Inputs: fname - input volume file corresponding to intput data for
%                 obtaining the voxel size
%         S_orig_data - the 3D input volume 
%         slices - slices that should be retained from interpolated volume
%         method - either 'fft' or interpolation method of 'interp3'
%         F - Interpolation factor (1 - no interpolation, >1 upsampling, ...)
% Output: Interpolated 3D volume
%
if ~sum(F_old - F_new)
    S = S_orig_data;
else
    S_orig_type = class(S_orig_data);
    if isempty(slices)
        slices = 1:size(S_orig_data, 3);
    end

    [X_old, Y_old, Z_old, X_new, Y_new, Z_new] = get_gridcoords(S_orig_data, F_old, F_new);
    Z_old_slices = Z_old(slices);
    M_z_new_slices = Z_new >= min(Z_old_slices) & Z_new <= max(Z_old_slices);
    Z_new_slices = Z_new(M_z_new_slices);

    switch method
        case 'fft'
            if ~sum(size(S_orig_data) - [length(X_old) length(Y_old) length(Z_old)])
                S_tmp = fftInterpolate( double(S_orig_data), ...
                                        [length(Y_new) length(X_new) length(Z_new)]);
                S = cast(S_tmp(:, :, M_z_new_slices), S_orig_type);
            else
                % currently does not support downsampling
                error('interp_series:input', 'Input data size missmatch!');
            end
        otherwise
            if ~sum(size(S_orig_data) - [length(X_old) length(Y_old) length(Z_old)])
                % prepare grids
                [Xo, Yo, Zo] = meshgrid(Y_old, X_old, Z_old_slices);
                [Xn, Yn, Zn] = meshgrid(Y_new, X_new, Z_new_slices);
                % prepare original data
                if size(S_orig_data, 3) == length(Z_old_slices)
                    S_tmp = double(S_orig_data);
                else
                    if size(S_orig_data, 3) > length(Z_old_slices)
                        S_tmp = double(S_orig_data(:, :, slices));
                    else
                        error('interp_series:input', 'Input data size missmatch!');
                    end
                end
            else
                if ~sum(size(S_orig_data) - [length(X_new) length(Y_new) length(Z_new_slices)])
                    % prepare grids
                    [Xo, Yo, Zo] = meshgrid(X_new, Y_new, Z_new_slices);    
                    [Xn, Yn, Zn] = meshgrid(X_old, Y_old, Z_old_slices);
                    % prepare original data
                    S_tmp = double(S_orig_data);
                else
                    error('interp_series:input', 'Input data size missmatch!');
                end
            end
            if isempty(method)
                method = 'linear';
            end
            S = cast(interp3(Xo, Yo, Zo, S_tmp, Xn, Yn, Zn, method), S_orig_type);
    end
end