function [SM_oli, CC] = morph_filter(SM_oli, S, S_roi, SM_ntis, SM_oli_hypo, F)

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
    S_out_1 = S_tmp(Roi.Maxrect.X, Roi.Maxrect.Y, Slices) .* 0;
    S_out_2 = S_tmp(Roi.Maxrect.X, Roi.Maxrect.Y, Slices) .* 0;
    for idx_slice = 1:N_slices
        S_in = S_tmp(Roi.Maxrect.X, Roi.Maxrect.Y, Slices(idx_slice));
        SM_in = S_in > 0;
        S_out_1(:, :, idx_slice) = nlfilter(S_in, [2 1], @get_contrast) .* cast(SM_in, class(S_in));
        S_out_2(:, :, idx_slice) = nlfilter(S_in, [1 2], @get_contrast) .* cast(SM_in, class(S_in));
    end
    S_out_1(isnan(S_out_1) | isinf(S_out_1) | S_out_1 == 0) = [];
    S_out_2(isnan(S_out_2) | isinf(S_out_2) | S_out_2 == 0) = [];
    Noise = quantile([S_out_1(:); S_out_2(:)], .05);
% end

fprintf('Before: %d ...', sum(SM_oli(:)));
L = conncomp_init(SM_oli, 3);
% S_gre_hdr = load_series(fname, 0);
% F = S_gre_hdr.hdr.dime.pixdim(2:4);
% S_gre = double(load_series(fname, slices));
% L_iso = interp_series(fname, L, slices, 'nearest', [1 1 1]);
CC = conncomp_list(L, L, S_roi, S, F, SM_oli_hypo);
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
    CC = CC(M);
end
% if CC(1).lab > -1
%     N_cc = length(CC);
% 	SM_oli = false(size(SM_oli));
%     idx_cc_new = 1;
%     for idx_cc = 1:N_cc
%         M_lab = Lab == CC(idx_cc).loc;
%         if CC(idx_cc).cont > Noise(M_lab)
%             SM_oli = SM_oli | L == CC(idx_cc).lab;
%             CC_new(idx_cc_new) = CC(idx_cc);
%             idx_cc_new = idx_cc_new + 1;
%         end
%     end
%     CC = CC_new;
% end
fprintf('After: %d.\n', sum(SM_oli(:)));


%-------------------------------------------------------------------------
function [cont] = get_contrast(S)
if sum(S(:)<=0)
    cont = NaN;
else
    Q = quantile(double(S(:)), [.25 .5 .75]);
    cont = (Q(3)-Q(1)) / Q(2);
end