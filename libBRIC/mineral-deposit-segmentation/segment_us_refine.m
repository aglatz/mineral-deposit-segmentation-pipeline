function [Ret] = segment_us_refine(Subject, RoiLabelTable, ReportName, ...
                                   InterpFactor, ThreshFactor, AdaptiveFlag, ...
                                   IntvarP)
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

Ret = validate(Subject, 'FE_roi_mask', 'FE_roi_mask', [], [], 0);
if Ret.Vol_ref > 0
    % Get optimal factor
    x0 = ThreshFactor(1);
    options = optimset('fminsearch');
    options = optimset(options, 'Display', 'iter');
    xopt = fminsearch(@(x) optifun(x, Subject, RoiLabelTable, ...
                      InterpFactor, AdaptiveFlag, IntvarP), x0, options);
                         
    % Get masks for best version
    Ret = segment_us_single(Subject, RoiLabelTable, ReportName, ...
                            InterpFactor, [xopt 0], AdaptiveFlag, IntvarP);
    Ret.alpha_opt = xopt;
else
	% Get masks for best version
    Ret = segment_us_single(Subject, RoiLabelTable, ReportName, ...
                            InterpFactor, ThreshFactor, AdaptiveFlag, IntvarP);
    Ret.alpha_opt = NaN;
end


%--------------------------------------------------------------------------
function [F] = optifun(x, Subject, RoiLabelTable, InterpFactor, AdaptiveFlag, IntvarP)
Ret = segment_us_single(Subject, RoiLabelTable, [], ...
                        InterpFactor, [x 0], AdaptiveFlag, IntvarP, ...
                        'SaveMaskFlag', false);
F = 1 - Ret.Jaccard;

