function [intvar_thr] = get_intvar_ntis(SM_ntis, S, intvar_q)
% Estimate T2*w intensity variance of normal-appearing tissue
% INPUTS: SM_ntis - Normal appearing tissue mask
%         S - T2*w intensities
%         intvar_q - Quantile for calulating the threshold
% OUTPUTS: intvar_thr - T2*w variance threshold associated with intvar_q
%

S_tmp = S .* cast(SM_ntis, class(S));
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
intvar_thr = quantile(IntVar, intvar_q);

%-------------------------------------------------------------------------
function [cont] = get_intvar_slicepatch(S)
if sum(S(:)<=0)
    cont = NaN;
else
    cont = var(double(S(:)));
end
