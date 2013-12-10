function [S_hypos, S_hypos_hypo, S_hypos_hyper, App] = ...
	cc_filter(S_gre, TE_gre, S_t1w, S_t2w, TE_t2w, S_roi, SM_hypos, ...
              S_ntis, I_thr, dR2s_thr, phypo_thr)
% This function filters the connected components of the preliminary
% basal ganglia T2*w hypointensity masks according to their median
% delta R2* value and appearance on T1w volumes. Connected components
% above dR2s_thr and below phypo_thr will be included in the final masks.
% INPUTS: S_gre - T2*w volume
%         TE_gre - Echo time of T2*w volume
%         S_t1w - Coregistered T1w volume
%         S_t2w - Experimental parameter
%         TE_t2w - Experimental parameter
%         S_roi - ROI mask set
%         SM_hypos - Preliminary basal ganglia T2*w hypointensity masks 
%         S_ntis - Normal-appearing tissue mask set
%         I_thr - T2*w and T1w thresholds for each ROI
%         dR2s_thr - delta R2* threshold
%         phypo_thr - T1w hypointensity threshold
% OUTPUTS: S_hypos - Final basal ganglia T2*w hypointensity mask
%          S_hypos_hypo - Regions of S_hypo with intensities below
%          5th-percentile
%          S_hypos_hyper - Regions of S_hypo with intensities below
%          95th-percentile

% Determine connected components (CCs)
L = conncomp_init(SM_hypos, 3);
[Lab, Loc] = conncomp_mask(L, S_roi, 0.5, S_gre);

fprintf('Before: %d; dR2s_thr=%0.2f phypo_thr=%0.2f\n', ...
        sum(SM_hypos(:)), dR2s_thr, phypo_thr);

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
    
    % Get R2* or R2' of features relative to 'normal-appearing' tissue
    % and check if the feature is Ok for inclusion in the final mask
    vol = sum(SM_cc(:));
    if isempty(S_t2w)
        phypo = get_phypo(SM_cc, S_t1w, I_thr(M, 4));
        phyper = get_phyper(SM_cc, S_t1w, I_thr(M, 5));
        dR2s = -1/TE_gre*log(median(S_gre(SM_cc))/median(S_gre(SM_ntis)));
        res = dR2s > dR2s_thr && vol > 1 && phypo < phypo_thr;
    else
        phypo = get_phypo(SM_cc, S_t2w, I_thr(M, 4));
        phyper = get_phyper(SM_cc, S_t2w, I_thr(M, 5));
        dR2s = -1/TE_gre*log(median(S_gre(SM_cc))/median(S_gre(SM_ntis))) - ...
               -1/TE_t2w*log(median(S_t2w(SM_cc))/median(S_t2w(SM_ntis)));
        res = dR2s > dR2s_thr && vol > 1 && phypo < phypo_thr && ...
              phyper < phypo_thr && ~(phypo > 0.05 && phyper > 0.05);
    end
    fprintf('Lab:%d loc:%d vol:%d phypo:%0.2f phyper:%0.2f dR2s:%0.2f', ...
        lab_idx, Loc(lab_idx), vol, phypo, phyper, dR2s);
	App(lab_idx, 1:(end-1)) = [sum(SM_cc(:)) Loc(lab_idx) vol phypo phyper dR2s];
    
    % Make decision ...
    if res
        App(lab_idx, end) = 1;
        % all
        S_hypos = S_hypos + cast(SM_cc, class(S_hypos)) .* Loc(lab_idx);
        % hypo hypo
        [~, SM_tmp] = get_phypo(SM_cc, S_t1w, I_thr(M, 2)); %*1.2);
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
