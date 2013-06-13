function [SM_oli] = morph_filter(SM_oli, S, S_roi, SM_ntis, SM_oli_hypo)

% Lab = unique(S_roi(:));
% Lab = Lab(2:end);
% N_lab = length(Lab);
% Noise = zeros(N_lab, 1);
% for idx_lab = 1:N_lab
%     S_tmp = S .* cast(S_roi == Lab(idx_lab) & SM_ntis, class(S));
S_tmp = S .* cast(logical(S_roi) & SM_ntis, class(S));
Roi = roi_init(S_tmp);
Slices = roi_nifti_sliceno(Roi, []);
N_slices = length(Slices);
S_out = S_tmp(Roi.Maxrect.X, Roi.Maxrect.Y, Slices) .* 0;
for idx_slice = 1:N_slices
    S_in = S_tmp(Roi.Maxrect.X, Roi.Maxrect.Y, Slices(idx_slice));
    SM_in = S_in > 0;
    S_out(:, :, idx_slice) = nlfilter(S_in, [2 2], @get_contrast) .* cast(SM_in, class(S_in));
end
S_out(isnan(S_out) | isinf(S_out) | S_out == 0) = [];
Noise = quantile(S_out, .25);
% end

fprintf('Before: %d ...', sum(SM_oli(:)));
L = conncomp_init(SM_oli, 3);
% S_gre_hdr = load_series(fname, 0);
% dx = S_gre_hdr.hdr.dime.pixdim(3);
% dy = S_gre_hdr.hdr.dime.pixdim(2);
% dz = S_gre_hdr.hdr.dime.pixdim(4);
% S_gre = double(load_series(fname, slices));
% L_iso = interp_series(fname, L, slices, 'nearest', [1 1 1]);
CC = conncomp_list(L, L, S_roi, S, [1 1 2], SM_oli_hypo);
if CC(1).lab > -1
    Noise'
    [CC.lab]
    [CC.cont]
    M = [CC.cont] > Noise; % & ~([CC.loc] ~= 13 & [CC.phypo] > 0.1 & [CC.vol] > 22);
    Lab = [CC(M).lab]';
    N_lab = length(Lab);
    SM_oli = false(size(SM_oli));
    for idx_lab = 1:N_lab
        SM_oli = SM_oli | L == Lab(idx_lab);
    end
end
fprintf('After: %d.\n', sum(SM_oli(:)));

%-------------------------------------------------------------------------
function [cont] = get_contrast(S)
if sum(S(:)<=0)
    cont = NaN;
else
    Q = quantile(double(S(:)), [.25 .5 .75]);
    cont = (Q(3)-Q(1)) / Q(2);
end