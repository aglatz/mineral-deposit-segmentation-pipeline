function [CC, SM, slices, L] = conncomp_stats_simple(path, name_con, use_lab, ...
                                                  name_roi, use_roi, use_iso, ...
                                                  name_weight)
if use_roi
    S_roi = load_series([path '/' name_roi], []);
    Roi = roi_init(S_roi);
	slices = roi_nifti_sliceno(Roi, []);
    % add one slice at the end so interpolation can work correctly
    % and we don't get the volume mismatch error below.
    slices_min = min(slices) - 1;
    if slices_min < 1
        slices_min = 1;
    end
    slices_max = max(slices) + 1;
    if slices_max > size(S_roi, 3)
        slices_max = size(S_roi, 3);
    end
    slices = (slices_min:slices_max)';
else
	slices = [];
end

% Get connected component labels
if use_lab
    % We get the labels from the concomp mask. That way we are consistent
    % with the result plots.
    L = load_series([path '/' name_con], slices);
    SM = logical(L);
else
	SM = logical(load_series([path '/' name_con], slices));
	L = conncomp_init(SM, 3);
end

% Get voxel size and isotropic voxel size
NII = load_series([path '/' name_roi], 0);
F = NII.hdr.dime.pixdim(2:4);
if use_iso
    F_iso = [min(F) min(F) min(F)]/2;
else
    F_iso = F;
end

% Get connected component labels of mask with isotropic voxels
if use_lab
    L_iso = load_series_interp([path '/' name_con], slices, ...
                               'nearest', F_iso);
else
    S_iso = load_series_interp([path '/' name_con], slices, ...
                               'nearest', F_iso);
%     S_iso = load_series([path '/' name_con], slices);
%     SM_tmp = isnan(S_iso);
%     S_iso(SM_tmp) = zeros(sum(SM_tmp(:)), 1);
	L_iso = conncomp_init(logical(S_iso), 3);
end

% Consistency check
if get_volume(logical(L), F) ~= get_volume(logical(L_iso), F_iso)
    error('Unequal voxel size!');
end

% Get Location labels and weights
S_roi = load_series([path '/' name_roi], slices);
S_weight = load_series([path '/' name_weight], slices);

% Get CC stats
CC = conncomp_list_simple(L, L_iso, S_roi, S_weight, F);

