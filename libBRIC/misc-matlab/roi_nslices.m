function [nslices] = roi_nslices(Roi)
% Return number of slices of ROI mask
nslices = size(Roi.Int, 3); % direct reference to Int because it's called
                            % from roi_mask()

