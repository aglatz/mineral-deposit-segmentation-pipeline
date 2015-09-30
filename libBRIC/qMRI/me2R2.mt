addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');
addpath('${SCRIPT_DIR}/LMFnlsq');

S_gre = load_series('${ARG_0}', []); % Multi echo gradient echo magnitude
N_echo = [ ${ARG_1} ];
S_gre = single(S_gre(:, :, :, 1:N_echo));
T_delta = ${ARG_3}-${ARG_2};
T = ${ARG_2}:T_delta:(${ARG_2}+T_delta*(N_echo-1)); % First echo and spacing in us
fprintf('Echoes: ');
for echo = 1:N_echo
	fprintf('%d ', T(echo));
end
fprintf('\n');

% Get brain mask
I_thr_input = str2num('${ARG_4}');
if isempty(I_thr_input)
	SM_brain = logical(load_series('${ARG_4}', [])); % Brain mask
	I_thr_input = [];
else
	SM_brain = S_gre(:, :, :, 1) > I_thr_input;
end

% Get cutoff threshold (no fitting if voxel intensity is below it)
S_tmp = S_gre(:, :, :, 1);
if isempty(I_thr_input)
	I_thr = quantile(S_tmp(SM_brain), .05);
else
	I_thr = I_thr_input;
end
fprintf('Lower threshold: %d\n', round(I_thr));

try
	% Estimate noise with dilated brain mask selecting
	% just the air (with noise)
	SM_brain_dil = logical(load_series('${ARG_5}', []));
	fprintf('Estimating noise...');
	Noise = zeros(N_echo, 1);
	for echo  = 1:N_echo
		S_tmp = S_gre(:, :, :, echo);
		Noise(echo) = iqr(S_tmp(SM_brain_dil))/(2*0.655);
		fprintf('%0.1f...', Noise(echo));
	end
	fprintf('done.\n');
catch
	Noise = ones(N_echo, 1);
end

% Recon
[R2s, R2s_sd, S0, S0_sd, R2s_csq, R2s_log] = ...
	recon_r2smap_lmf(S_gre, Noise, T'./1e6, I_thr, 1);

% Save results
[path, name, ext] = fileparts('${ARG_6}');
if isempty(path)
	path = '.';
end
if strcmp(ext, '.gz')
	[tmp, name, ext] = fileparts(name);
	ext = [ext '.gz'];
end
save_series('${ARG_0}', fullfile(path, name), single(R2s), []);
save_series('${ARG_0}', fullfile(path, [name '_std']), single(R2s_sd), []);
save_series('${ARG_0}', fullfile(path, 'T1w'), single(S0), []);
save_series('${ARG_0}', fullfile(path, ['T1w_std']), single(S0_sd), []);
save_series('${ARG_0}', fullfile(path, [name '_csq']), single(R2s_csq), []);
save_series('${ARG_0}', fullfile(path, [name '_log']), int16(R2s_log), []);
S_tmp = R2s_csq(SM_brain);
S_tmp(isnan(S_tmp)) = [];
SM = R2s_csq > quantile(S_tmp, .05) & R2s_csq < quantile(S_tmp, .95);
save_series('${ARG_0}', fullfile(path, [name '_sos_90']), SM, []);
SM = R2s_csq > quantile(S_tmp, .1) & R2s_csq < quantile(S_tmp, .9);
save_series('${ARG_0}', fullfile(path, [name '_sos_80']), SM, []);

