function [Ret] = segment_us(S_gre, S_t1w, S_voi, S_ref, Lab, ...
                            I_ntis_means_est, I_thr, adaptive_flag, varargin)
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
%         I_thr - If not empty it is used for the segmentation and
%                 circumvents the calculation of the adaptive threshold
%         adaptive_flag - 'true' indicates that the adaptive outlier
%                         detection method should be used
% OPTIONAL INPUTS: Out_name - File name of the ps file where the plots
%                             summary the performance statistics are saved.
% RETURN STRUCTURE ELEMENTS:
%          S_hypos - The basal ganglia hypointensity masks
%          S_nontis - The outlier masks
%          S_ntis - The masks of normal-appearing tissues
%          I_ntis_means - GRE and T1W robust mean intensities of the 
%                         normal-appearing tissue.
%          I_thr - The derived segmentation thresholds
%

% ROI labels - Threshold I_gre_min_min is derived from the first label
N_lab = length(Lab);
% Colors for the plot: black - background, green - ok, red - not ok, blue - missing
Col = [[0, 0, 0]; hsv(3)];
% Return structure
Ret = struct;
% Masks that selects the hypointense outlier intensities
Ret.S_hypos = zeros(size(S_voi), class(S_voi));
% Mask that selects tissue intensities not belonging to a specific label
Ret.S_nontis = zeros(size(S_voi), class(S_voi));
% Mask selecting normal tissue
Ret.S_ntis = zeros(size(S_voi), class(S_voi));
% Robust normal tissue means
Ret.I_ntis_means = zeros(N_lab, 3);
% Segmentation thresholds
Ret.I_thr = NaN(N_lab, 6);
% Fit
Ret.Fit = NaN(N_lab, 1);
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
            case {'S_ref_ca'}
                S_ref_ca = varargin{idx_vain+1};    
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
    if exist('S_ref_ca', 'var') && ~isempty(S_ref_ca)
        SM_ref = SM_ref & logical(S_ref_ca);
    end
    
    % Figure for scatter and chi2 plots
    H = figure;

    % Plot all tissue intensities
    subplot(2, 1, 1);
    psize = 20;
    scatter(S_gre(SM_voi), S_t1w(SM_voi), psize, Col(1, :));
    title(sprintf('ID=%d', Lab(idx_lab)));
    hold on;
    
	% Get means and covariance matrix of approximately normally 
    % distributed intensities (the mode of the GMM, which is closest
    % to I_ntis_mean_est or the biggest cluster)
    if ~isempty(I_ntis_means_est)
        I_ntis_mean_est = I_ntis_means_est(idx_lab, :);
        scatter(I_ntis_mean_est(1), I_ntis_mean_est(2), 20, 'r', 'filled');
    else
        I_ntis_mean_est = [];
    end
    [I_ntis_means, C_ntis, SM_mode, GM_2] = ...
        get_antissue(S_gre, S_t1w, SM_voi, I_ntis_mean_est, 1);
    Ret.I_ntis_means(idx_lab, :) = [I_ntis_means, Lab(idx_lab)];
    
    % Plot approximately normally distributed intensities
    if ~isempty(GM_2)
%         % Discrimination line
%         a0=log(GM_2.PComponents(1)/GM_2.PComponents(2)) - 1/2 * ...
%             ((GM_2.mu(1, :)+GM_2.mu(2, :)) / GM_2.Sigma) * (GM_2.mu(1, :)-GM_2.mu(2, :))';
%         a1 = GM_2.Sigma \ (GM_2.mu(1, :)-GM_2.mu(2, :))';
%         P=[-a1(1)/a1(2) -a0/a1(2)];
%         X = min(S_gre(SM_mode)):10:max(S_gre(SM_mode));
%         Y = polyval(P, X);
%         plot(X, Y, 'r');
%         text(round(max(X)), round(max(Y)), sprintf('y = %0.2f x  + %0.2f', P(1), P(2)));
        
        % Selected cluster
        scatter(S_gre(SM_mode), S_t1w(SM_mode), psize, 'c');
        scatter(I_ntis_means(1), I_ntis_means(2), 20, 'm', 'filled');
    end
    
    % Get normal tissue, candidate outliers and thresholds
    [SM_oli, SM_ntis, RDs] = ...
        get_normal_outliers(S_gre, S_t1w, SM_voi, I_ntis_means, C_ntis, adaptive_flag);
    Ret.S_nontis = Ret.S_nontis + cast(SM_oli, class(Ret.S_nontis)) .* Lab(idx_lab);
    Ret.S_ntis = Ret.S_ntis + cast(SM_ntis, class(Ret.S_ntis)) .* Lab(idx_lab);
    
    % Thresholding (Derive threshold from first ROI if not given)
    if idx_lab == 1
        [SM_hypos, I_thr] = ...
            thresh_filter(S_gre, S_t1w, SM_oli, I_ntis_means, C_ntis, ...
                          RDs, P_thr, []);
    else
        [SM_hypos, I_thr] = ...
            thresh_filter(S_gre, S_t1w, SM_oli, I_ntis_means, C_ntis, ...
                          RDs, P_thr, I_thr);      
    end
    Ret.I_thr(idx_lab, :) = [I_thr, Lab(idx_lab)];

    % Save hypointensity masks
    Ret.S_hypos = Ret.S_hypos + cast(SM_hypos, class(Ret.S_hypos)) .* Lab(idx_lab);
    
    subplot(2, 1, 1);
	% Plot missclassification error
    scatter(S_gre(SM_hypos), S_t1w(SM_hypos), psize, Col(2, :), 'filled');
    scatter(S_gre(SM_ref & SM_hypos), S_t1w(SM_ref & SM_hypos), psize, Col(3, :), 'filled');
	scatter(S_gre(SM_ref & ~SM_hypos), S_t1w(SM_ref & ~SM_hypos), psize, Col(4, :), 'filled');
    xlabel('\bf T2*W in arb. units');
    ylabel('\bf T1W in arb. units');
    % Show thresholds
    ellipsplot_mod(I_ntis_means, C_ntis, [], 'b', RDs(1));
    ellipsplot_mod(I_ntis_means, C_ntis, [], 'g', RDs(2));
    vline(I_thr(1), 'k', '');
    hline(I_thr(2), 'k--', '');
    hline(I_thr(3), 'k--', '');
    hline(I_thr(4), 'k:', '');
    hline(I_thr(5), 'k:', '');

    % Plot mahalanobis against chi2
    subplot(2, 1, 2);
    plot_mahalanobis(S_gre, S_t1w, SM_voi, I_ntis_means, C_ntis, RDs);
    
    % save figure
	save_ps_figure(Out_name, H);
    
    % Distributions
    H = figure;
    subplot(3, 1, 1);
    Ret.Fit(idx_lab) = plot_rddist(S_gre, S_t1w, SM_voi, SM_ntis, I_ntis_means, C_ntis, RDs);
    
    subplot(3, 1, 2);
    plot_intdist(S_gre, SM_voi, SM_ntis, I_thr(1), I_ntis_means(1), 'T2*w');
    
    subplot(3, 1, 3);
	plot_intdist(S_t1w, SM_voi, SM_ntis, I_thr(2:5), I_ntis_means(2), 'T1w');
    
	% save figure
	save_ps_figure(Out_name, H);
end


%% Plot robust distance distribution
function [P_fit] = plot_rddist(S_gre, S_t1w, SM_voi, SM_ntis, I_ntis_mean, C_ntis, RDs)
% 
Mat = [S_gre(SM_ntis) S_t1w(SM_ntis)];
%Mat = mvnrnd(I_ntis_mean, C_ntis, 1000); % Test

MD_ntis_oli_ref = (mahalanobis(Mat, I_ntis_mean, 'cov', C_ntis));
MD_ntis_oli_ref(MD_ntis_oli_ref >= RDs(2)) = [];
[Y1, X] = plot_hist(MD_ntis_oli_ref, [], [], 'b', 0); % dummy - just to get the values

% Sample truncated chi2 distribution
% https://stat.ethz.ch/pipermail/r-help/2012-February/302876.html
N_samp = length(MD_ntis_oli_ref);
B = chi2inv(rand(N_samp, 1)*chi2cdf(RDs(2), 2), 2);

% Plot MDs and Chi2s
Y2 = hist(B, X);
hold off;
plot_bar2(X, Y1, Y2, 'RD', 'Chi2');
hold on;
xlabel('\bf Robust distances in arb. units');
ylabel('\bf Occurrence');

% Test for equal population
[H, P_fit] = kstest2(MD_ntis_oli_ref, B);
%[P_fit, H] = ranksum(MD_ntis_oli_ref, B);
title(sprintf('N_total=%d, N_norm=%d, kstest2={H=%d, P=%0.3f}', sum(SM_voi(:)), sum(SM_ntis(:)), H, P_fit));


%% Plot intensity distribution
function plot_intdist(S, SM_voi, SM_antis, Thr, DistMean, ti)
[Y1, X] = plot_hist(S(SM_voi), [], [], 'g', 0);
Y2 = hist(S(SM_antis), X);
hold off;
plot_bar2(X, Y1, Y2, 'VOI', 'Nor');
hold on;
xlabel(['\bf ' ti ' in arb. units']);
ylabel('\bf Occurrence');
title(sprintf('N=%d', sum(SM_antis(:))));
if ~isempty(Thr)
    N_thr = length(Thr);
    for idx = 1:N_thr
        vline(Thr(idx), 'g', '');
    end
end
vline(DistMean, 'b', '');
% Bowley skewness
Q = quantile(S(SM_voi), [.25 .5 .75]);
Skew = (Q(1)-2*Q(2)+Q(3))/(Q(3)-Q(1));
title(sprintf('IQR=%0.3f skew=%0.3f', (Q(3)-Q(1))/Q(2), Skew));


%% Plot mahalanobis distance against sqrt of chi2 distribution
function plot_mahalanobis(S_gre, S_t1w, SM_voi, I_ntis_mean, C_ntis, RDs)
% QQ plot
Mat = [S_gre(SM_voi) S_t1w(SM_voi)];
MD_ntis_oli_ref = (mahalanobis(Mat, I_ntis_mean, 'cov', C_ntis));
% [y, x] = chiqqplot_mod(MD_ntis_oli_ref, 2, 'MCDCOV');
N_samp = length(MD_ntis_oli_ref);
if N_samp < 1000
    N_samp = 1000;
end
h = qqplot((chi2rnd(2, N_samp, 1)), MD_ntis_oli_ref);
x = get(h(1), 'XData');
y = get(h(1), 'YData');
set(h(1), 'LineStyle', '-');
set(h(1), 'Marker', 'none');
xlabel('\bf Quantiles of the \chi^2 distribution (df=2)');
ylabel('\bf Robust distances');
hold on;

% Display fixed threshold
diff = (get(h(1), 'YData') - (RDs(1))).^2;
y_idx = find(min(diff) == diff, 1);
scatter(x(y_idx), y(y_idx), 'b', 'filled');

% Display refined threshold
diff = (y - (RDs(2))).^2;
y_idx = find(min(diff) == diff, 1);
scatter(x(y_idx), y(y_idx), 'g', 'filled');


%% Morphological filtering to remove outliers from segmentation, imaging and vessels 
function [SM_hypo, I_thr] = ...
    thresh_filter(S_gre, S_t1w, SM_oli, I_ntis_mean, C_ntis, RDs, P_thr, I_thr)
% Derive threshold if not given
if isempty(I_thr)
    % Calculate GRE threshold with refined RD
    I_gre_min = -sqrt(C_ntis(1,1)*RDs(2))+I_ntis_mean(1);
    I_thr(1) = polyval(P_thr, I_gre_min);
    fprintf('Thr=%0.3f(diff=%0.3f)...', I_thr, abs(I_thr(1)-I_ntis_mean(1)));
end

% Calculate T1W thresholds
I_thr(2) = -sqrt(C_ntis(2,2))*1.645+I_ntis_mean(2); % 90% range
I_thr(3) = +sqrt(C_ntis(2,2))*1.645+I_ntis_mean(2);
I_thr(4) = -sqrt(C_ntis(2,2)*RDs(2))+I_ntis_mean(2); % Robust range
I_thr(5) = +sqrt(C_ntis(2,2)*RDs(2))+I_ntis_mean(2);

% Calculate delta delta R2 dash
% dR2s = -1/15e-3*log(I_thr(1)/I_ntis_mean(1));
% dR2 = -1/102.96e-3*log(I_thr(4)/I_ntis_mean(2));
% ddR2d = dR2s - dR2;
% fprintf('ddR2s=%0.3f...', ddR2d);

% Segment multifocal T2*w hypointensities
Mat = [S_gre(SM_oli) S_t1w(SM_oli)];
SM_hypo = SM_oli;
SM_hypo(SM_oli) = Mat(:, 1) < I_thr(1); % Apply GRE threshold
fprintf('T2*w hypos: %d(diffmin=%0.3f).\n', sum(SM_hypo(:)), abs(max(S_gre(SM_hypo))-I_ntis_mean(1)));



