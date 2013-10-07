function [ddR2d] = get_ddR2smap(S_gre, TE_gre, S_t1w, TE_t1w, S_roi, I_thr)
Labs = unique(S_roi(:));
Labs = Labs(2:end); % Exclude background
N_labs = length(Labs);
ddR2d = NaN(size(S_gre));
for idx_lab = 1:N_labs
    SM_roi = S_roi == Labs(idx_lab);
    M = I_thr(:, end) == Labs(idx_lab);
    ddR2d(SM_roi) = -1/TE_gre*log(S_gre(SM_roi)./I_thr(M, 1)) - ...
                (-1/TE_t1w*log(S_t1w(SM_roi)./I_thr(M, 2)));
end