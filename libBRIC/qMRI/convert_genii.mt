addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

% Input params
fname_in = '${ARG_0}';
N_te = ${ARG_1};
te = ${ARG_2}; % in us
te2 = ${ARG_3}; % in us
if ${ARG_4} > 0 && ${ARG_5} > 0
	Res = [${ARG_4} ${ARG_5}];
else
	Res = [];
end
fname_out = '${ARG_6}';

% Calculate TEs
if N_te > 1
    T = (te:(te2-te):te+(te2-te)*(N_te-1))';
else
    T = te;
end

% Load multi-echo data
[S_abs, S_pha] = read_genii(fname_in, T, ${ARG_7}, Res);

% Save multi-echo data
save_series([fname_in '_' num2str(T(1))], [fname_out '_mag'], S_abs, []);
save_series([fname_in '_' num2str(T(1))], [fname_out '_pha'], S_pha, []);
save_series([fname_in '_' num2str(T(1))], [fname_out '_rea'], real(S_abs.*exp(sqrt(-1).*S_pha)), []);
save_series([fname_in '_' num2str(T(1))], [fname_out '_ima'], imag(S_abs.*exp(sqrt(-1).*S_pha)), []);

% Save k-space data
save_series_kspace(S_abs, S_pha, fname_out, T);

% Save ln of mag data
S_tmp = log(S_abs);
SM = isnan(S_tmp) | isinf(S_tmp);
S_tmp(SM) = 0;
save_series([fname_in '_' num2str(T(1))], [fname_out '_mag_ln'], S_tmp, []);

