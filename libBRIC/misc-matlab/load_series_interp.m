function S = load_series_interp(fname, slices, method, F_new)
% Loads 3D volume and interpolates it to the given voxel size.
% Inputs: fname - input volume
%         slices - slices that should be loaded from input volume
%         method - either 'fft' or interpolation method of 'interp3'
%         F_new - New voxel size
%
if isempty(F_new)
    S = load_series(fname, slices);
else
    S_orig_data = load_series(fname, []);
    NII = load_series(fname, 0);
    F_old = NII.hdr.dime.pixdim(2:4);
    S = interp_series(S_orig_data, slices, method, F_old, F_new);
end

       
