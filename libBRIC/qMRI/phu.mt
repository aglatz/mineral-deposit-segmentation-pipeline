addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/LIBRA');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

% Input data
fname_in = '${ARG_0}';
N_te = ${ARG_1};
te = double(${ARG_2})./1e6; % in us
te2 = double(${ARG_3})./1e6; % in us
fname_out = '${ARG_7}';

% Calculate TEs
if N_te > 1
    T = (te:(te2-te):te+(te2-te)*(N_te-1))';
else
    T = te;
end

% Load input data
S_abs = double(load_series([fname_in '_mag'], []));
S_abs = S_abs(:, :, :, (1:N_te)); %*2-1);
S_pha = double(load_series([fname_in '_pha'], []));
S_pha = S_pha(:, :, :, (1:N_te)); %*2-1);

% Get brain mask
I_thr_input = str2num('${ARG_5}');
if isempty(I_thr_input)
	SM_brain = logical(load_series('${ARG_5}', [])); % Brain mask
	I_thr_input = [];
else
	SM_brain = S_abs(:, :, :, 1) > I_thr_input;
end

try
	% Estimate noise
	SM_noise = logical(load_series('${ARG_6}', []));
	Noise = zeros(N_te, 1);
	for idx = 1:N_te
		S_tmp = S_abs(:, :, :, idx);
		Noise(idx) = iqr(S_tmp(SM_noise))/(2*0.655);
	end
catch
	Noise = ones(N_te, 1);
end
Noise

% Remove low frequency variations
NII = load_series([fname_in '_mag'], 0);
F = NII.hdr.dime.pixdim(2:4);
for idx=1:N_te
	S_tmp = S_abs(:, :, :, idx) .* exp(sqrt(-1)*S_pha(:, :, :, idx));
	N = size(S_tmp);
	S_tmp_f = gauss3filter(padarray(S_tmp, round(N/2)), ${ARG_4}, F);
	S_tmp = S_tmp ./ S_tmp_f(round(N(1)/2)+1:round(N(1)/2)+N(1), ...
                        	 round(N(2)/2)+1:round(N(2)/2)+N(2), ...
                        	 round(N(3)/2)+1:round(N(3)/2)+N(3));
	S_abs(:, :, :, idx) = abs(S_tmp) .* cast(SM_brain, class(S_tmp));
	S_pha(:, :, :, idx) = angle(S_tmp) .* cast(SM_brain, class(S_tmp));
end
save_series([fname_in '_mag'], [fname_out '_filt'], single(S_pha), []);

% Reconstruct unwrapped phase by assuming a linear phase evolution: phi = k*TE + d
K_est = zeros(size(S_pha(:, :, :, 1)));
D_est = zeros(size(S_pha(:, :, :, 1)));
SoS_est = NaN(size(S_pha(:, :, :, 1)));
for z=1:size(K_est, 3)
    fprintf('%d ...', z);
    for y=1:size(K_est, 2)
        for x=1:size(K_est, 1)
			if SM_brain(x, y, z)
            	I = unwrap(reshape(S_pha(x, y, z, :), N_te, 1));
            	if ~sum(isnan(I))
					A = [T(:) ones(N_te, 1)];
					b = I(:);
					w = 1./Noise(:).^2;
					[P, ~, MSE] = lscov(A, b, w);
					K_est(x, y, z) = P(1);
					D_est(x, y, z) = P(2);
					SoS_est(x, y, z) = MSE;
            	end
			end
        end
    end
end
fprintf('done.\n');

% Save result
save_series([fname_in '_mag'], [fname_out '_phuf0'], ...
			single(K_est./(2*pi)), []);
save_series([fname_in '_mag'], [fname_out '_sos'], SoS_est, []);
S_tmp = SoS_est(SM_brain);
S_tmp(isnan(S_tmp)) = [];
SM = SoS_est > quantile(S_tmp, .05) & SoS_est < quantile(S_tmp, .95);
save_series([fname_in '_mag'], [fname_out '_sos_90'], SM, []);
SM = SoS_est > quantile(S_tmp, .1) & SoS_est < quantile(S_tmp, .9);
save_series([fname_in '_mag'], [fname_out '_sos_80'], SM, []);

