function [ResMat, N_res] = get_cc_avg_stats_all(CC, Loc, Q, Field, cc_type)
% Convenience function that calls 'get_cc_avg_stats' for each field
% and location given as input.
% Inputs: CC - list of structures
%         Loc - location labels of connected components that should be
%               analyzed
%         Q - quantiles (e.g. [25 50 75])
%         Field - fields to summarize
%         cc_type - see 'get_cc_avg_stats'
% Outputs: ResMat - Quantiles for each field and location combination
%          N_res - Count
%
N_Loc = length(Loc);
N_Field = length(Field);
ResMat = zeros(N_Field, N_Loc, length(Q));
N_res = zeros(N_Loc, 1);
for idx_field = 1:N_Field
    for idx_loc = 1:N_Loc
        [N, Val]  = get_cc_avg_stats(CC, Loc{idx_loc}, Field{idx_field}, cc_type, '');
        N_res(idx_loc) = sum(N);
        ResMat(idx_field, idx_loc, :) = quantile(cell2mat(Val), Q);
    end
end
