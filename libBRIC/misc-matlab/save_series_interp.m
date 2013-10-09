function [] = save_series_interp(fname, new_name, S_orig_data, slices, method, F_new)
% Interpolates the given volume to the given voxel size and saves it.
% Inputs: fname - template volume
%         new_name - new volume name
%         S_orig_data - Data to save
%         slices - slices that were loaded from original volume or empty
%         method - interpolation method of 'interp3'
%         F_new - voxel size of interpolated volume
%
if isempty(F_new)
    save_series(fname, new_name, S_orig_data, slices);
else
    NII = load_series(fname, 0);
    F_old = NII.hdr.dime.pixdim(2:4);
    S = interp_series(S_orig_data, slices, method, F_new, F_old);
    save_series(fname, new_name, S, slices);
end

