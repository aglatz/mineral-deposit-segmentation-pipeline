function [S_out, S_oli_all, S_ntis, I_ntis_means, I_gre_thr, I_gre_thr_ref] =...
    segment_us(S_gre, S_t1w, S_voi, S_ref, Lab, I_ntis_means_est, I_gre_thr, ...
               varargin)
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
% RETURNS: S_out - The basal ganglia hypointensity masks
%          S_oli_all - The outlier masks
%          S_ntis - The masks of normal-appearing tissues
%          I_ntis_means - GRE and T1W robust mean intensities of the 
%                         normal-appearing tissue.
%          I_gre_thr - The derived GRE segmentation thresholds
%          I_gre_thr_ref - The GRE reference thresholds
%

% ROI labels - Threshold I_gre_min_min is derived from the first label
N_lab = length(Lab);
% Colors for the plot: black - background, green - ok, red - not ok, blue - missing
Col = [[0, 0, 0]; hsv(3)];
% Masks that selects the hypointense outlier intensities
S_out = zeros(size(S_voi), class(S_voi));
% Mask that selects tissue intensities not belonging to a specific label
S_oli_all = zeros(size(S_voi), class(S_voi));
% Mask selecting normal tissue
S_ntis = zeros(size(S_voi), class(S_voi));
% Robust normal tissue means
I_ntis_means = zeros(N_lab, 2);
% Chi2 for given area under Chi2 distribution with 2 deg freedom
Chi2_thr = [];
% Linear polynomial for adjusting the threshold
P_thr = [1 0];
% Minimum contrast of a connected component
cont_min = -0.001;
% Name of file were plots get saved to
Out_name = 'segment_us_plots';
% Process optional args
N_vain = length(varargin);
for idx_vain = 1:N_vain
    arg_in = varargin{idx_vain};
    if ischar(arg_in) || isscalar(arg_in)
        switch arg_in
            case {'Chi2_thr'}
                Chi2_thr = varargin{idx_vain+1};
            case {'Out_name'}
                Out_name = varargin{idx_vain+1};
            case {'P_thr'}
                P_thr = varargin{idx_vain+1};
        end
    end
end

% Reference intensity percentiles
Pctls = [.7 .75 .8 .85 .9];
I_gre_thr_ref = NaN(1, length(Pctls));

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
    H = create_ps_figure;  %figure;

    % Plot all tissue intensities
    subplot(2, 1, 1);
    scatter(S_gre(SM_voi), S_t1w(SM_voi), 10, Col(1, :));
    hold on;
    
	% Get approximately normally distributed intensities
    [SM_antis, GM_2] = get_antissue(S_gre, S_t1w, SM_voi, I_ntis_means_est(idx_lab, :), 2);
    
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
    [SM_oli, SM_ntis, I_gre_min, I_t1w_min, I_ntis_means(idx_lab, :), C_ntis] = ...
        get_normal_outliers(S_gre, S_t1w, SM_antis, SM_voi, Chi2_thr);
    S_oli_all = S_oli_all + cast(SM_oli, class(S_oli_all)) .* Lab(idx_lab);
    S_ntis = S_ntis + cast(SM_ntis, class(S_ntis)) .* Lab(idx_lab);
    
    % Filter out outliers from segmentation, imaging and vessels.
    if isempty(I_gre_thr)
        % Derive threshold from first ROI
        if Lab(idx_lab) == Lab(1)
            % Derive Threshold from GP
            [SM_oli_filt, I_gre_thr] = ...
                morph_filter(S_gre, S_t1w, SM_oli, I_gre_min, I_t1w_min, [], P_thr, 1, SM_voi, cont_min);
            
            % Save reference GP intensity percentiles
            if sum(SM_ref(:)) > 0
                I_gre_thr_ref = quantile(S_gre(SM_ref), Pctls);
            end
        else
            [SM_oli_filt] = ...
                morph_filter(S_gre, S_t1w, SM_oli, I_gre_min, I_t1w_min, I_gre_thr, [], 1, SM_voi, cont_min);
        end
    else
        [SM_oli_filt] = ...
                morph_filter(S_gre, S_t1w, SM_oli, I_gre_min, I_t1w_min, I_gre_thr, [], 0, SM_voi, cont_min);
    end
    
    % Plot missclassification error
    subplot(2, 1, 1);
    scatter(S_gre(SM_oli_filt), S_t1w(SM_oli_filt), 10, Col(2, :));
    scatter(S_gre(SM_ref & SM_oli_filt), S_t1w(SM_ref & SM_oli_filt), 10, Col(3, :));
	scatter(S_gre(SM_ref & ~SM_oli_filt), S_t1w(SM_ref & ~SM_oli_filt), 10, Col(4, :));
    xlabel('\bf T2*W in arb. units');
    ylabel('\bf T1W in arb. units');
    h(1) = line_with_label(I_gre_thr, 'thr=', 'g', max(S_t1w(SM_voi)), [], 1, 'v');
    I_tmp = quantile(S_gre(SM_ref), .95);
    h(2) = line_with_label(I_tmp, 'thr=', '--b', max(S_t1w(SM_voi)), [], 1, 'v');
    legend(h, {'Estimated', 'Reference'}, 'Location', 'best');

    % Plot mahalanobis against chi2
    subplot(2, 1, 2);
    plot_mahalanobis(S_gre, S_t1w, SM_ntis, SM_oli_filt, SM_ref, I_ntis_means(idx_lab, :), C_ntis);

    % Save result
    S_out = S_out + cast(SM_oli_filt, class(S_out)) .* Lab(idx_lab);
    
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
function [SM_oli, I_gre_thr] = ...
    morph_filter(S_gre, S_t1w, SM_oli, I_gre_min, I_t1w_min, I_gre_thr, P_thr, adj_thr, SM_voi, cont_min)

Mat = [S_gre(SM_oli) S_t1w(SM_oli)];

% Derive threshold if not given
if isempty(I_gre_thr)
    I_gre_thr = polyval(P_thr, I_gre_min);
    fprintf('%0.3f(%0.3f)...', I_gre_thr, I_gre_min);
end

% Adjust thresholds if necessary - must be outside tolerance ellipse
if adj_thr && I_gre_thr > I_gre_min
    I_gre_thr = I_gre_min;
end

% apply T1W&GRE thresholds (soft)
SM_tmp_red = SM_oli;
SM_tmp_red(SM_oli) = Mat(:, 1) < I_gre_thr & Mat(:, 2) > I_t1w_min;
% apply GRE threshold (hard)
SM_oli(SM_oli) = Mat(:, 1) < I_gre_thr;

fprintf('Before: %d ...', sum(SM_oli(:)));
% Fill holes
L = conncomp_init(SM_oli, 3);
Lab = unique(L)';
Lab = Lab(2:end);
Lab_red = unique(L(SM_tmp_red))';
N_lab = length(Lab);
SM_oli = false(size(SM_tmp_red));
for idx_lab = 1:N_lab
    if sum(ismember(Lab(idx_lab), Lab_red))
        SM_oli = SM_oli | L == Lab(idx_lab);
    end
end
% Filter out low contrast CC
L = conncomp_init(SM_oli, 3);
CC = conncomp_list(L, L, SM_voi, S_gre, [1 1 1]);
if CC(1).lab > -1
    M = [CC.cont] < cont_min;
    Lab = [CC(M).lab]';
    N_lab = length(Lab);
    SM_oli = false(size(SM_oli));
    for idx_lab = 1:N_lab
        SM_oli = SM_oli | L == Lab(idx_lab);
    end
end
fprintf('After: %d.\n', sum(SM_oli(:)));


%% Split normal and outlier intensities
function [SM_oli, SM_ntis, I_gre_min, I_t1w_min, I_ntis_mean, C_ntis] = ...
    get_normal_outliers(S_gre, S_t1w, SM_antis, SM_voi, Chi2_thr)

Mat = [S_gre(SM_antis) S_t1w(SM_antis)];

% Mean and covariance of approx. normal distr. tissue
[~, I_ntis_mean, I_ntis_sd, ~, ~, C_ntis] = pcomp_find(Mat);
if isempty(Chi2_thr)
    Chi2_thr = get_chi2_thr(2, size(Mat, 1));
end

% Split normal tissue and outliers
Mat = [S_gre(SM_voi) S_t1w(SM_voi)];
Chi2 = get_chi_right_tail(Mat', I_ntis_mean, C_ntis);
M = Chi2(:) < Chi2_thr;
SM_ntis = SM_voi;
SM_ntis(SM_voi) = M;
SM_oli = SM_voi;
SM_oli(SM_voi) = ~M;

% Thresholds
[I_gre_min, I_t1w_min] = ellipsplot_mod(I_ntis_mean, C_ntis, Mat, 'r', Chi2_thr);
pcomp_plot(I_ntis_mean, I_ntis_sd, 'b');


%%
function [chi2_thr] = get_chi2_thr(p, n)
chi2_thr = chi2inv(1-(0.24-0.003*p)/sqrt(n), p);
