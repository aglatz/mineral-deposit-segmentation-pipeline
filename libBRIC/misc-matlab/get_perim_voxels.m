function [SM_perim, A] = get_perim_voxels(SM, interal)
% Returns a mask that selects the perimeter voxels of the ROI mask.
% INPUTS SM - ROI mask
%        internal - flag indicating if the perimeter mask should
%                   select voxels that surround the ROI mask (and
%                   are not part of the mask) or the voxels that 
%                   are the edge of the ROI mask.
% OUTPUTS SM_perim - the perimeter mask
%         A - number of voxels of the perimeter mask
%
e1 = [[0 0 0];[0 1 0]; [0 0 0]];
e2 = [[0 1 0];[1 0 1]; [0 1 0]];
SE = cat(3, e1, e2);
SE = cat(3, SE, e1);
I_tmp = ones(size(SM))*6 - convn(SM, SE, 'same');
if interal
    SM_perim = SM & convn(SM, SE, 'same') < 6;
else
    SM_perim = ~SM & convn(SM, SE, 'same') > 0;
end
A = get_volume(I_tmp .* cast(SM_perim, 'double'));
