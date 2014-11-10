%% Split normal and outlier intensities
function [SM_oli, SM_ntis, RDs] = ...
    get_normal_outliers(S_gre, S_t1w, SM_voi, I_ntis_mean, C_ntis, adaptive_flag)

Mat = [S_gre(SM_voi) S_t1w(SM_voi)];

% Get robust distances of selected mode
RD = mahalanobis(Mat, I_ntis_mean, 'cov', C_ntis);
delta = chi2inv(0.975, size(Mat, 2)); % Only simulation values for 97.5%
if adaptive_flag
    % Find outliers according to Filzmoser, Reimann, Garrett (2003).
    [Gn, u] = ecdf(RD);
    G = chi2cdf(u, size(Mat, 2));
    dG = G - Gn;
    pn = max(dG(u >= delta & dG >= 0));
    an = 0;
    if ~isempty(pn)
        pcrit = (0.24 - 0.003 * size(Mat, 2))/sqrt(size(Mat, 1));
        if pn > pcrit
            an = pn;
        end    
    end
    RD_cutoff = u(find(Gn > (1-an), 1));
    if isempty(RD_cutoff)
        RD_cutoff = max(RD);
    end
else
    RD_cutoff = delta;
end
RDs = [delta RD_cutoff];

% Thresholding
M = RD <= RD_cutoff;
SM_ntis = SM_voi;
SM_ntis(SM_voi) = M;
SM_oli = SM_voi;
SM_oli(SM_voi) = ~M;

