function [] = save_series_interp(fname, new_name, S_orig_data, slices, method, F)
% Interpolates the given volume to the given voxel size and saves it.
% Inputs: fname - template volume
%         new_name - new volume name
%         S_orig_data - Data to save
%         slices - slices that were loaded from original volume or empty
%         method - interpolation method of 'interp3'
%         F - Interpolation factor (1 - no interpolation, >1 upsampling, ...)
% Outputs: S - Interpolated volume
%
S = interp_series(fname, S_orig_data, slices, method, F);
save_series(fname, new_name, S, slices);
