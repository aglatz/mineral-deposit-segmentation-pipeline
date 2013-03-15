function [sliceno] = roi_nifti_sliceno(Roi, sliceidx)
% Takes structure returned by roi_init() and returns the
% slice index of all slices covered by the ROI.
%
if isempty(sliceidx)
    sliceidx = 1:roi_nslices(Roi);
end
sliceno = Roi.Num(sliceidx, 1);

