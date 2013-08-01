function [SM_oli_out, CC, IntVar] = intvar_filter(SM_oli_in, S, S_roi, SM_ntis, ...
                                                  SM_oli_hypo, F, intvar_p)
% This function estimates the T2*w intensity variance of normal-appearing
% tissue selected by SM_ntis and calculates the intensity variablity of
% the connected components of the SM_oli_in mask representing focal T2*w
% hypointensities. Connected components are retained if their intensity
% variability is higher than the estimated one of normal tissue.
% INPUTS: SM_oli_in - Preliminary focal T2*w hypointensity mask
%         S - T2*w volume
%         S_roi - ROI mask
%         SM_ntis - Mask that selects only normal-appearing tissue
%         SM_oli_hypo - Mask selecting voxels that appear hypointense on
%                       T1w volumes
%         F - Voxel size
%         intvar_p - Cummulative probability value that defines the
%                    T2*w intensity variability threshold
%
% OUTPUTS: SM_oli_out - Filtered focal T2*w hypointensity mask
%          CC - Structure that lists the statistics of connected components
%          IntVar - Vector containing normal-appearing intensity
%                   variabilities
%
                                              
% Estimate T2*w intensity variance of normal-appearing tissue
S_tmp = S .* cast(logical(S_roi) & SM_ntis, class(S));
Roi = roi_init(S_tmp);
Slices = roi_nifti_sliceno(Roi, []);
N_slices = length(Slices);
S_out = S_tmp(Roi.Maxrect.X, Roi.Maxrect.Y, Slices) .* 0;
for idx_slice = 1:N_slices
    S_in = S_tmp(Roi.Maxrect.X, Roi.Maxrect.Y, Slices(idx_slice));
    SM_in = S_in > 0;
    S_out(:, :, idx_slice) = nlfilter(S_in, [2 2], @get_intvar_slicepatch) ...
                             .* cast(SM_in, class(S_in));
end
S_out(isnan(S_out) | isinf(S_out) | S_out == 0) = [];
IntVar = S_out(:);
PD = fitdist(IntVar, 'lognormal');
intvar_thr = PD.icdf(intvar_p);

% Get T2*w intensity variability of individual T2*w hypointensities and
% exclude the ones below the noise threshold
fprintf('Before: %d ...', sum(SM_oli_in(:)));
L = conncomp_init(SM_oli_in, 3);
CC = conncomp_list(L, L, S_roi, S, F, SM_oli_hypo);
if CC(1).lab > -1
    intvar_thr'
    [CC.lab]
    [CC.intvar]
	M = [CC.intvar] > intvar_thr;
    Lab = [CC(M).lab]';
    N_lab = length(Lab);
    SM_oli_out = false(size(SM_oli_in));
    for idx_lab = 1:N_lab
        SM_oli_out = SM_oli_out | L == Lab(idx_lab);
    end
    CC = CC(M);
else
    SM_oli_out = SM_oli_in;
end
fprintf('After: %d.\n', sum(SM_oli_out(:)));


%-------------------------------------------------------------------------
function [cont] = get_intvar_slicepatch(S)
if sum(S(:)<=0)
    cont = NaN;
else
    cont = var(double(S(:)));
end