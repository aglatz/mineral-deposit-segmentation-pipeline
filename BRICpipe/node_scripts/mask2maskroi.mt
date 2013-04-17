addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');
	
SM = logical(load_series('${ARG_0}', []));
L = conncomp_init(SM, ${ARG_1});
S_roi = load_series('${ARG_2}', []);
S_weight = load_series('${ARG_3}', []);
[Lab, Loc] = conncomp_mask(L, S_roi, ${ARG_4}, S_weight);
M = Loc > -1;
Lab = Lab(M);
N_lab = length(Lab);
Loc = Loc(M);
S_filter = zeros(size(S_roi), class(S_roi));
for idx = 1:N_lab
    SM = L == Lab(idx);
    S_filter = S_filter + Loc(idx) * cast(SM, class(S_filter));
end
save_series('${ARG_2}', '${ARG_5}', S_filter, []);

