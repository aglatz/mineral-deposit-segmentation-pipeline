function [Ret] = segment_us(S_gre, S_t1w, S_voi, S_ref, Lab, ...
                            I_ntis_means_est, I_gre_thr, varargin)
% Unsupervised segmentation of hypointense outlier signal intensities on GRE,
% which are likely iron deposits or calcifications. Requires GRE and T1-weighted
% volumes as input. T1-weighted volumes are used to separate hypointensities
% that are predominatly caused by iron deposits from hypointensities that are
% predominantly caused by calcifications. Processes all VOI regions with labels
% given in vector 'Lab'. An adaptive threshold is derived from the VOI region
% with label 'Lab(1)' and used for the remaining ROI regions. Usually the first
% VOI region is the globus pallidus and the other regions are the caudate, 
% putamen, and internal capsule. This function also requires the mean 
% tissue intensities of the normal-appearing tissue of all structures. The
% adaptive thresholding is circumvented if 'I_gre_thr' is not empty. Plots
% are generated that can be used to analyze the performance of the
% segmentation method. Those plots also contain info about the voxels selected
% by the reference mask 'S_ref' if not empty.
% INPUTS: S_gre - GRE volume
%         S_t1w - Co-registered T1W volume
%         S_voi - Combined VOI masks, where each mask has a unique label
%         S_ref - Used for deriving performance statistics if not empty
%         Lab - VOI mask labels corresponding to the regions that should
%               be analyzed. The adaptive threshold is derived from the
%               first region and used for the segmentation of the other
%               regions.
%         I_ntis_means_est - Mean GRE and T1W intensities of the normal-
%                            appearing tissue of all regions.
%         I_gre_thr - If not empty it is used for the segmentation and
%                     circumvents the calculation of the adaptive threshold
% OPTIONAL INPUTS: Out_name - File name of the ps file where the plots
%                             summary the performance statistics are saved.
% RETURN STRUCTURE ELEMENTS:
%          S_out - The basal ganglia hypointensity masks
%          S_nontis - The outlier masks
%          S_ntis - The masks of normal-appearing tissues
%          I_ntis_means - GRE and T1W robust mean intensities of the 
%                         normal-appearing tissue.
%          I_gre_thr - The derived GRE segmentation thresholds
%

% ROI labels - Threshold I_gre_min_min is derived from the first label
N_lab = length(Lab);
% Colors for the plot: black - background, green - ok, red - not ok, blue - missing
Col = [[0, 0, 0]; hsv(3)];
% Return structure
Ret = struct;
% Masks that selects the hypointense outlier intensities
Ret.S_out = zeros(size(S_voi), class(S_voi));
Ret.S_out_hypo = zeros(size(S_voi), class(S_voi));
Ret.S_out_hyper = zeros(size(S_voi), class(S_voi));
% Mask that selects tissue intensities not belonging to a specific label
Ret.S_nontis = zeros(size(S_voi), class(S_voi));
% Mask selecting normal tissue
Ret.S_ntis = zeros(size(S_voi), class(S_voi));
% Robust normal tissue means
Ret.I_ntis_means = zeros(N_lab, 2);
% Linear polynomial for adjusting the threshold
P_thr = [1 0];
% Name of file were plots get saved to
Out_name = 'segment_us_plots';
% Process optional args
N_vain = length(varargin);
for idx_vain = 1:N_vain
    arg_in = varargin{idx_vain};
    if ischar(arg_in) || isscalar(arg_in)
        switch arg_in
            case {'Out_name'}
                Out_name = varargin{idx_vain+1};
            case {'P_thr'}
                P_thr = varargin{idx_vain+1};
        end
    end
end

% --- Main loop ---
for idx_lab = 1:N_lab
	fprintf('Processing label: %d...', Lab(idx_lab));
	SM_voi = S_voi == Lab(idx_lab);
    if ~isempty(S_ref)
        SM_ref = S_ref == Lab(idx_lab);
    else
        SM_ref = false(size(SM_voi));
    end
    
    % Figure for scatter and chi2 plots
    H = figure; %create_ps_figure;  %figure;

    % Plot all tissue intensities
    subplot(2, 1, 1);
    scatter(S_gre(SM_voi), S_t1w(SM_voi), 10, Col(1, :));
    hold on;
    
	% Get approximately normally distributed intensities
    if isempty(I_ntis_means_est)
        I_ntis_mean_est = [];
    else
        I_ntis_mean_est = I_ntis_means_est(idx_lab, :);
    end
%     [SM_antis, GM_2] = get_antissue(S_gre, S_t1w, SM_voi, I_ntis_mean_est, 2);
    [SM_antis, GM_2] = get_antissue(S_gre, S_t1w, SM_voi, I_ntis_mean_est, 1);
    
    % Plot approximately normally distributed intensities
    if ~isempty(GM_2)
%         % Discrimination line
%         a0=log(GM_2.PComponents(1)/GM_2.PComponents(2)) - 1/2 * ...
%             ((GM_2.mu(1, :)+GM_2.mu(2, :)) / GM_2.Sigma) * (GM_2.mu(1, :)-GM_2.mu(2, :))';
%         a1 = GM_2.Sigma \ (GM_2.mu(1, :)-GM_2.mu(2, :))';
%         P=[-a1(1)/a1(2) -a0/a1(2)];
%         X = min(S_gre(SM_antis)):10:max(S_gre(SM_antis));
%         Y = polyval(P, X);
%         plot(X, Y, 'r');
%         text(round(max(X)), round(max(Y)), sprintf('y = %0.2f x  + %0.2f', P(1), P(2)));
        
        % Selected cluster
        scatter(S_gre(SM_antis), S_t1w(SM_antis), 10, 'c');
    end
    
    % Get normal tissue, candidate outliers and thresholds
    [SM_oli, SM_ntis, I_gre_min, I_t1w_min, I_t1w_max, Ret.I_ntis_means(idx_lab, :), C_ntis] = ...
        get_normal_outliers(S_gre, S_t1w, SM_antis, SM_voi);
    Ret.S_nontis = Ret.S_nontis + cast(SM_oli, class(Ret.S_nontis)) .* Lab(idx_lab);
    Ret.S_ntis = Ret.S_ntis + cast(SM_ntis, class(Ret.S_ntis)) .* Lab(idx_lab);
    
    % Filter out outliers from segmentation, imaging and vessels.
    if isempty(I_gre_thr)
        % Derive threshold from first ROI
        if Lab(idx_lab) == Lab(1)
            % Derive Threshold from GP
            [SM_oli_filt, Ret.I_gre_thr, SM_oli_filt_t1whypo, SM_oli_filt_t1whyper] = ...
                thresh_filter(S_gre, S_t1w, SM_oli, I_gre_min, I_t1w_min, I_t1w_max, [], P_thr, 1);
        else
            [SM_oli_filt, tmp, SM_oli_filt_t1whypo, SM_oli_filt_t1whyper] = ...
                thresh_filter(S_gre, S_t1w, SM_oli, I_gre_min, I_t1w_min, I_t1w_max, Ret.I_gre_thr, [], 1);
        end
    else
        [SM_oli_filt, Ret.I_gre_thr, SM_oli_filt_t1whypo, SM_oli_filt_t1whyper] = ...
                thresh_filter(S_gre, S_t1w, SM_oli, I_gre_min, I_t1w_min, I_t1w_max, I_gre_thr, [], 0);
    end
    
    % Plot missclassification error
    subplot(2, 1, 1);
    scatter(S_gre(SM_oli_filt), S_t1w(SM_oli_filt), 10, Col(2, :));
    scatter(S_gre(SM_ref & SM_oli_filt), S_t1w(SM_ref & SM_oli_filt), 10, Col(3, :));
	scatter(S_gre(SM_ref & ~SM_oli_filt), S_t1w(SM_ref & ~SM_oli_filt), 10, Col(4, :));
    xlabel('\bf T2*W in arb. units');
    ylabel('\bf T1W in arb. units');
    I_tmp = quantile(S_gre(SM_oli_filt), .5);
    h(1) = line_with_label(I_tmp, 'thr=', 'g', max(S_t1w(SM_voi)), [], 1, 'v');
    I_tmp = quantile(S_gre(SM_ref), .5);
    h(2) = line_with_label(I_tmp, 'thr=', '--b', max(S_t1w(SM_voi)), [], 1, 'v');
    legend(h, {'Estimated median', 'Reference median'}, 'Location', 'best');

    % Plot mahalanobis against chi2
    subplot(2, 1, 2);
    plot_mahalanobis(S_gre, S_t1w, SM_ntis, SM_oli_filt, SM_ref, Ret.I_ntis_means(idx_lab, :), C_ntis);

    % Save result
    Ret.S_out = Ret.S_out + cast(SM_oli_filt, class(Ret.S_out)) .* Lab(idx_lab);
    Ret.S_out_hypo = Ret.S_out_hypo + cast(SM_oli_filt_t1whypo, class(Ret.S_out_hypo)) .* Lab(idx_lab);
    Ret.S_out_hyper = Ret.S_out_hyper + cast(SM_oli_filt_t1whyper, class(Ret.S_out_hyper)) .* Lab(idx_lab);
    
    % save figure
	save_ps_figure(Out_name, H);
end


%% Plot mahalanobis distance against sqrt of chi2 distribution
function [MD_oli_min, MD_ntis_oli, MD_min_ref] = ...
    plot_mahalanobis(S_gre, S_t1w, SM_ntis, SM_oli, SM_ref, I_ntis_mean, C_ntis)

% QQ plot
SM_ntis_oli_ref = SM_ntis | SM_oli | SM_ref;
Mat = [S_gre(SM_ntis_oli_ref) S_t1w(SM_ntis_oli_ref)];
MD_ntis_oli_ref = sqrt(mahalanobis(Mat, I_ntis_mean, 'cov', C_ntis));
[y, x] = chiqqplot_mod(MD_ntis_oli_ref, 2, 'MCDCOV');
hold on;

% Plot chosen threshold
Leg = {};
idx = 1;
MD_oli_min = min(MD_ntis_oli_ref(SM_oli(SM_ntis_oli_ref)));
if ~isempty(MD_oli_min)
    idx_y = find(y == MD_oli_min, 1);
    h(1) = scatter(x(idx_y), MD_oli_min, 'g', 'filled');
    Leg{idx} = 'Estimated';
    idx = idx+1;
end
% Plot reference threshold
if ~isempty(SM_ref)
    MD_min_ref = min(MD_ntis_oli_ref(SM_ref(SM_ntis_oli_ref)));
    if ~isempty(MD_min_ref)
        idx_y = find(y == MD_min_ref, 1);
        h(2) = scatter(x(idx_y), MD_min_ref, 'b', 'filled');
        Leg{idx} = 'Reference';
        % idx = idx+1;
    end
else
    MD_min_ref = [];
end
if false && ~isempty(Leg)
	legend(h, Leg, 'Location', 'best');
end

% Plot adjusted outlyingness
SM_ntis_oli = SM_ntis | SM_oli;
MD_ntis_oli = MD_ntis_oli_ref(SM_ntis_oli(SM_ntis_oli_ref));


%% Morphological filtering to remove outliers from segmentation, imaging and vessels 
function [SM_oli, I_gre_thr, SM_oli_t1whypo, SM_oli_t1whyper] = ...
    thresh_filter(S_gre, S_t1w, SM_oli, I_gre_min, I_t1w_min, I_t1w_max, ...
                  I_gre_thr, P_thr, adj_thr)
% Derive threshold if not given
if isempty(I_gre_thr)
    I_gre_thr = polyval(P_thr, I_gre_min);
    fprintf('%0.3f(%0.3f)...', I_gre_thr, I_gre_min);
end

% Adjust thresholds if necessary - must be outside tolerance ellipse
if adj_thr && I_gre_thr > I_gre_min
    I_gre_thr = I_gre_min;
end

% Morphological reconstruction
Mat = [S_gre(SM_oli) S_t1w(SM_oli)];
SM_oli_marker = SM_oli; % Marker volume: apply T1W&GRE thresholds (soft)
SM_oli_marker(SM_oli) = Mat(:, 1) < I_gre_thr ...
                        & Mat(:, 2) > I_t1w_min & Mat(:, 2) < I_t1w_max;
SM_oli(SM_oli) = Mat(:, 1) < I_gre_thr; % Mask volume: apply GRE threshold (hard)
fprintf('Before morph. recon.: %d ...', sum(SM_oli(:)));
SM_oli = imreconstruct(SM_oli_marker, SM_oli, 6);
fprintf('After: %d.\n', sum(SM_oli(:)));

% Save regions that are hypo- or hyperintense on T1w
Mat = [S_gre(SM_oli) S_t1w(SM_oli)];
SM_oli_t1whypo = SM_oli;
SM_oli_t1whypo(SM_oli) = Mat(:, 2) <= I_t1w_min;
SM_oli_t1whyper = SM_oli;
SM_oli_t1whyper(SM_oli) = Mat(:, 2) >= I_t1w_max;



%% Split normal and outlier intensities
function [SM_oli, SM_ntis, I_gre_min, I_t1w_min, I_t1w_max, I_ntis_mean, C_ntis] = ...
    get_normal_outliers(S_gre, S_t1w, SM_antis, SM_voi)

Mat = [S_gre(SM_antis) S_t1w(SM_antis)];

% Mean and covariance of approx. normal distr. tissue
[~, I_ntis_mean, I_ntis_sd, ~, ~, C_ntis] = pcomp_find(Mat);

% Find outliers according to Filzmoser, Reimann, Garrett (2003).
RD = mahalanobis(Mat, I_ntis_mean, 'cov', C_ntis);
[Gn, u] = ecdf(RD);
G = chi2cdf(u, size(Mat, 2));
delta = chi2inv(0.975, size(Mat, 2)); % Only simulation values for 97.5%
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
% RD_cutoff = delta;

% Thresholding
M = RD <= RD_cutoff;
SM_ntis = SM_voi;
SM_ntis(SM_voi) = M;
SM_oli = SM_voi;
SM_oli(SM_voi) = ~M;

% Show thresholds
ellipsplot_mod(I_ntis_mean, C_ntis, Mat, 'b', delta);
ellipsplot_mod(I_ntis_mean, C_ntis, Mat, 'r', RD_cutoff);
pcomp_plot(I_ntis_mean, I_ntis_sd, 'b');

% Calculate exact thresholds
I_gre_min = -sqrt(C_ntis(1,1)*RD_cutoff)+I_ntis_mean(1);
I_t1w_min = -sqrt(C_ntis(2,2)*RD_cutoff)+I_ntis_mean(2);
I_t1w_max = +sqrt(C_ntis(2,2)*RD_cutoff)+I_ntis_mean(2);

