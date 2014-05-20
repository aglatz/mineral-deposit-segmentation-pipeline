function [xopt] = segment_us_refine(Subjects, RoiLabelTable, N_cpus, AdaptiveFlag)
% Function that finds the optimal 'ThreshFactor' (see segment_us())
% by adjusting the factor to maximise the Jaccard index between
% generated and reference mask.
% INPUTS: ... - see segment_us_single()
% OUTPUTS:  ... - same as segment_us_single()
%           alpha_opt - the optimal slope of the ThreshFactor; the intercept
%                 is held at 0 and actually not needed. alpha_opt=NaN means
%                 that no reference mask was available and alpha was not
%                 optimized.
%

% Get optimal factor
intstd_thr_0 = 1;
options = optimset('fminsearch');
options = optimset(options, 'Display', 'iter');
xopt = fminsearch(@(x) optifun(Subjects, RoiLabelTable, N_cpus, ...
                  AdaptiveFlag, x), intstd_thr_0, options);


%--------------------------------------------------------------------------
function [F] = optifun(Subjects, RoiLabelTable, N_cpus, AdaptiveFlag, intstd_thr)
[Ret] = segment_us_mp(Subjects, RoiLabelTable, N_cpus, ...
                      'ThreshFactor', [1 0], 'AdaptiveFlag', AdaptiveFlag, ...
                      'SaveMaskFlag', false, 'N_gre', 1, ...
                      'CNR_thr', 0, 'phypo_thr', 0.1, 'intvar_thr', intstd_thr);
F = 1 - quantile([Ret.Jaccard], .5);
fd = fopen('/tmp/opt.txt', 'a');
fprintf(fd, '%s thr=%0.3f 1-J=%0.3f\n', datestr(now), intstd_thr, F);
fclose(fd);

