function [Mask] = roi_mask(Roi, sliceidx)
% Returns Roi slices that are non-zero in ROI mask
if isempty(sliceidx)
    sliceidx = 1:roi_nslices(Roi);
end
Mask = logical(Roi.Int(:, :, sliceidx));


