function [I_antis_mean, C_antis, SM_antis, GM_2] = ...
            get_antissue(S_gre, S_t1w, SM_roi, I_ntis_mean_est, N_comp_max)
% Applies the logical mask 'SM_roi' to the volumes 'S_gre' and 'S_t1w', 
% to obtain a bivariate T2*- and T1-weighted signal intensity distribution.
% Then a unimodal, multivariate Gaussian distribution and N=2...N_comp_max
% Gaussian mixture models (shared covariance matrix) are fitted to the data.
% The best fitting mixture model is chosen according to the minimal BIC
% (Fraley, JASA, 2002). The multivariate model is chosen over the univariate
% model if AIC is below 0 (Freeman, Behavior Research Methods, 2012). If 
% the multivariate model is chosen then the signal intensity cluster 
% that is closes to 'I_ntis_mean_est' is classified as normal appearing tissue.
% INPUTS: S_gre - T2*-weighted signal intensities
%         S_t1w - T1-weighted signal intensities
%         SM_roi - Logical ROI mask
%         I_ntis_mean_est - estimated normal tissue mean intensities
%         N_comp_max - max. number of Gaussian mixture components
% RETURNS: I_antis_mean - Robust mean T2*w/T1w tissue intensities
%          C_antis - Robust covariance matrix
%          SM_antis - Reduced logical ROI mask that just selects
%                     approximately normally distributed tissue intensities
%          GM_2 - Matlab structure of the fitted Gaussian mixture model
%          I_ntis_means - Normal-appearing tissue signal intensity means
%                         that were calculated with the 'SM_roi' mask.
%
Mat = [S_gre(SM_roi) S_t1w(SM_roi)];

% Find out whether distribution is multimodal
GM_1 = gmdistribution.fit(Mat, 1);
AIC_1 = GM_1.AIC;
N_comp = 1;
N_cnt = 3;
AICd = NaN;
cnt = 1;
while GM_1.Converged && cnt <= N_cnt && isnan(AICd)
    BIC_old = NaN;
    GM_2 = [];
    while N_comp < N_comp_max
        GM = gmdistribution.fit(Mat, N_comp+1, 'SharedCov', true);
        Idx = cluster(GM, Mat);
        Idx_u = unique(Idx)';
        N_Idx_u = length(Idx_u);
        N_el = zeros(N_Idx_u, 1);
        for idx = 1:N_Idx_u
            N_el(idx) = sum(Idx==Idx_u(idx));
        end
        if ~GM.Converged || (~isnan(BIC_old) && GM.BIC > BIC_old) || min(N_el) < 5
            break;
        end
        N_comp = N_comp+1;
        BIC_old = GM.BIC;
        GM_2 = GM;
    end
    if ~isempty(GM_2) && GM_2.Converged
        AIC_2 = GM_2.AIC;
        AICd = (AIC_1 - AIC_2)/max([AIC_1 AIC_2]);
    else
        AICd = NaN;
    end
    cnt = cnt + 1;
end
fprintf('#comp:%d(%d)...', N_comp, cnt-1);
if ~isnan(AICd) && AICd > 0
    Idx = cluster(GM_2, Mat);
    Idx_u = unique(Idx)';
    N_Idx_u = length(Idx_u);
    
    % Find biggest cluster
    N = zeros(N_Idx_u, 1);
    for idx = 1:N_Idx_u
        M = Idx == Idx_u(idx);
        N(idx) = sum(M);
    end
    M = N == max(N);

    % Generate mask
    SM_antis = SM_roi;
	SM_antis(SM_roi) = Idx == Idx_u(M);
	fprintf('#%d...', sum(SM_antis(:)));
else
    SM_antis = SM_roi;
    GM_2 = [];
end
% Mean and covariance of approx. normal distr. tissue
Mat = [S_gre(SM_antis) S_t1w(SM_antis)];
[~, I_antis_mean, ~, ~, ~, C_antis] = pcomp_find(Mat);
