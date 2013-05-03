addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');
addpath('${SCRIPT_DIR}/LMFnlsq');

S_gre = load_series('${ARG_0}', []); % Multi echo gradient echo magnitude
N_echo = [ ${ARG_1} ];
S_gre = single(S_gre(:, :, :, 1:N_echo));
T_delta = ${ARG_3}-${ARG_2};
T = ${ARG_2}:T_delta:(${ARG_2}+T_delta*N_echo); % First echo and spacing in us
fprintf('Echoes: ');
for echo = 1:N_echo
	fprintf('%d ', T(echo));
end
fprintf('\n');
I_thr_input = str2num('${ARG_4}');
if isempty(I_thr_input)
	SM_brain = logical(load_series('${ARG_4}', [])); % Brain mask
	I_thr_input = [];
end
SM_brain_dil = logical(load_series('${ARG_5}', [])); % Dilated brain mask
													 % selecting just the noise

% Estimate noise
fprintf('Estimating noise...');
Noise = zeros(N_echo, 1);
for echo  = 1:N_echo
	Noise_tmp = zeros(size(SM_brain_dil, 3), 1);
	for slice = 1:size(SM_brain_dil, 3)
		S_tmp = S_gre(:, :, slice, echo);
		% Noise extimation from magnitude requires a correction factor
		Noise_tmp(slice) = iqr(S_tmp(SM_brain_dil(:, :, slice)))/(2*0.655);
	end
	Noise(echo) = quantile(Noise_tmp, .5);
	fprintf('%0.1f...', Noise(echo));
end
fprintf('done.\n');
    
% Get lower cutoff
S_tmp = S_gre(:, :, :, 1);
if isempty(I_thr_input)
	I_thr = quantile(S_tmp(SM_brain), .10);
else
	I_thr = I_thr_input;
end
fprintf('Lower threshold: %d\n', round(I_thr));

% Recon
[R2s, R2s_sd, Tmp, Tmp, R2s_csq, R2s_log] = recon_r2smap_lmf(S_gre, Noise, T'./1e6, I_thr, 1);

% Save results
[path, name, ext] = fileparts('${ARG_6}');
if isempty(path)
	path = '.';
end
if strcmp(ext, '.gz')
	[tmp, name, ext] = fileparts(name);
	ext = [ext '.gz'];
end
save_series('${ARG_5}', [path '/' name], single(R2s), []);
save_series('${ARG_5}', [path '/' name '_std'], single(R2s_sd), []);
save_series('${ARG_5}', [path '/' name '_csq'], single(R2s_csq), []);
save_series('${ARG_5}', [path '/' name '_log'], int16(R2s_log), []);

