function S = load_series_interp(fname, slices, method, F)
% Loads 3D volume and interpolates it to the given voxel size.
% Inputs: fname - input volume
%         slices - slices that should be loaded from input volume
%         method - either 'fft' or interpolation method of 'interp3'
%         F - Interpolation factor (1 - no interpolation, >1 upsampling, ...)
%
S_orig_data = load_series(fname, []);
S = interp_series(fname, S_orig_data, slices, method, F);


       
