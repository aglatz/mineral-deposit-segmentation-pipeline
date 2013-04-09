function [Size] = roi_XYsize(Roi)
% Return size of ROI mask
Size = size(Roi.Int(:, :, 1));
