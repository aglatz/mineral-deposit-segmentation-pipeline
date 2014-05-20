function [S_hypos, S_hypos_hypo, S_hypos_hyper, App] = ...
	cc_filter(S_gre, S_t1w, N_gre, S_roi, SM_hypos, ...
              S_ntis, I_thr, CNR_thr, phypo_thr, intvar_thr)
% This function filters the connected components of the preliminary
% basal ganglia T2*w hypointensity masks according to their median
% delta R2* value and appearance on T1w volumes. Connected components
% above dR2s_thr and below phypo_thr will be included in the final masks.
% INPUTS: S_gre - T2*w volume
%         S_t1w - Coregistered T1w volume
%         N_gre - Estimated noise of T2*w volume
%         S_roi - ROI mask set
%         SM_hypos - Preliminary basal ganglia T2*w hypointensity masks 
%         S_ntis - Normal-appearing tissue mask set
%         I_thr - T2*w and T1w thresholds for each ROI from outlier
%                 detection method
%         CNR_thr - contrast-to-noise ratio threshold
%         phypo_thr - T1w hypointensity threshold
%         intvar_thr - T2*w intensity variance threshold
% OUTPUTS: S_hypos - Final basal ganglia T2*w hypointensity mask
%          S_hypos_hypo - Regions of S_hypo that appear hypointense on T1w
%          S_hypos_hyper - Regions of S_hypo that appear hyperintense on T1w
%          App - Apperance statistics (internal)

% Determine connected components (CCs)
L = conncomp_init(SM_hypos, 3);
[Lab, Loc] = conncomp_mask(L, S_roi, 0.5, S_gre);

% T2*w intensity variance of normal appearing tissue
S_gre_std = stdfilt(double(S_gre));

fprintf('Before: %d; CNR_thr=%0.2f phypo_thr=%0.2f instd_thr=%0.2f\n', ...
        sum(SM_hypos(:)), CNR_thr, phypo_thr, intvar_thr);

% Iterate through all CCs
N_lab = length(Lab);
S_hypos = zeros(size(S_roi), class(S_roi));
S_hypos_hypo = zeros(size(S_roi), class(S_roi));
S_hypos_hyper = zeros(size(S_roi), class(S_roi));
App = NaN(N_lab, 7);
for lab_idx = 1:N_lab
    % 3D mask of feature
	SM_cc = L == Lab(lab_idx);
    % 3D mask of 'normal-appearing' tissue
    SM_ntis = S_ntis == Loc(lab_idx);
    % Mask for getting the right T1w thresholds for discriminating outliers...
    M = I_thr(:, end) == Loc(lab_idx);

    % Get standardized T2*w intensity variance of connected-components
    intvar = std(double(S_gre(SM_cc)))/median(S_gre_std(SM_ntis));
    
    % Check if the connected-component is OK for inclusion in the final mask
    vol = sum(SM_cc(:));
    % M = I_thr(:, end) == 11; 
    phypo = get_phypo(SM_cc, S_t1w, I_thr(M, 4));
    % M = I_thr(:, end) == 14;
    phyper = get_phyper(SM_cc, S_t1w, I_thr(M, 5));
    CNR = abs(median(S_gre(SM_cc))-median(S_gre(SM_ntis)))/N_gre; %CNR
    res = CNR > CNR_thr && vol > 1 && phypo < phypo_thr && intvar > intvar_thr;
    fprintf('Lab:%d loc:%d vol:%d phypo:%0.2f phyper:%0.2f CNR:%0.2f intvar:%0.2f', ...
            lab_idx, Loc(lab_idx), vol, phypo, phyper, CNR, intvar);
	App(lab_idx, 1:(end-1)) = [sum(SM_cc(:)) Loc(lab_idx) vol phypo phyper CNR];
    
    % Make decision ...
    if res
        App(lab_idx, end) = 1;
        % all
        S_hypos = S_hypos + cast(SM_cc, class(S_hypos)) .* Loc(lab_idx);
        % hypo hypo
        [~, SM_tmp] = get_phypo(SM_cc, S_t1w, I_thr(M, 2));
        S_hypos_hypo = S_hypos_hypo + cast(SM_tmp, class(S_hypos_hypo)) .* Loc(lab_idx);
        % hypo hyper
        [~, SM_tmp] = get_phyper(SM_cc, S_t1w, I_thr(M, 3));
        S_hypos_hyper = S_hypos_hyper + cast(SM_tmp, class(S_hypos_hyper)) .* Loc(lab_idx);
        fprintf('*\n');
    else
        App(lab_idx, end) = 0;
        fprintf('\n');
    end
end

fprintf('After: %d\n', sum(logical(S_hypos(:))));


%-------------------------------------------------------------------------
function [phypo, SM_hypo] = get_phypo(SM, S, I_thr)
SM_hypo = false(size(SM));
SM_hypo(SM) = S(SM) < I_thr;
phypo = sum(SM_hypo(:))/sum(SM(:));

%-------------------------------------------------------------------------
function [phyper, SM_hyper] = get_phyper(SM, S, I_thr)
SM_hyper = false(size(SM));
SM_hyper(SM) = S(SM) > I_thr;
phyper = sum(SM_hyper(:))/sum(SM(:));

