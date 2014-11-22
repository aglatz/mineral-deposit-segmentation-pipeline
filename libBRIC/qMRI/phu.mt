addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/LIBRA');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');
addpath('${SCRIPT_DIR}/STI_Suite_v1.42/Core_Functions')

% ARG_0: Basename in
% ARG_1: Number of echoes
% ARG_2: TE1 in us
% ARG_3: TE2 in us
% ARG_4: Echo increment
% ARG_5: FWHM in mm
% ARG_6: == 0 -> low-pass, != 0 -> high-pass
% ARG_7: ICV mask
% ARG_8: Basename out 

% Input data
fname_in = '${ARG_0}';
skip = ${ARG_4};
N_te = ${ARG_1}/skip;
te = double(${ARG_2})./1e6; % in us
te2 = double(${ARG_3})./1e6; % in us
fname_out = '${ARG_8}';
fprintf('hu: %s\n', fname_out);
% Calculate TEs
if N_te > 1
    idx_te = skip:skip:N_te*skip;
    T = (te*skip:(te2-te)*skip:te*skip+(te2-te)*skip*(N_te-1))';
else
    idx_te = 1;
    T = te;
end
T

% Load input data
S_abs_in = double(load_series([fname_in '_mag'], []));
S_pha_in = double(load_series([fname_in '_pha'], []));
N_xy = [size(S_abs_in, 1), size(S_abs_in, 2), 1];
if mod(size(S_abs_in, 3), 2) == 1
	S_abs_in = cat(3, S_abs_in(:, :, :, idx_te), zeros([N_xy N_te]));
	S_pha_in = cat(3, S_pha_in(:, :, :, idx_te), zeros([N_xy N_te]));
	expanded = 1;
else
	S_abs_in = S_abs_in(:, :, :, idx_te);
	S_pha_in = S_pha_in(:, :, :, idx_te);
	expanded = 0;
end 

% Get brain mask
I_thr_input = str2num('${ARG_7}');
if isempty(I_thr_input)
	SM_brain = logical(load_series('${ARG_7}', [])); % Brain mask
	I_thr_input = [];
else
	SM_brain = S_abs_in(:, :, :, 1) > I_thr_input;
end
if expanded
	SM_brain = cat(3, SM_brain, false(N_xy));
end

% Remove low frequency variations
NII = load_series([fname_in '_mag'], 0);
F = NII.hdr.dime.pixdim(2:4);
N_pad = ceil(F*10);
sigma = ${ARG_5}/2.3548; % ARG_5 is FWHM
fprintf('gauss3filter: FWHM:%0.3fmm Sigma:%0.3fmm\n', ${ARG_5}, sigma);
S_phu = zeros(size(S_pha_in));
S_phu_f = zeros(size(S_pha_in));
for idx=1:N_te
	S_phu(:, :, :, idx) = LaplacianPhaseUnwrap(S_pha_in(:,:,:,idx), F, N_pad);
	S_tmp_pad = padarray(S_phu(:, :, :, idx), N_pad);
	S_tmp_pad_f = gauss3filter(S_tmp_pad, sigma, F);
	if ${ARG_6}
		S_tmp_pad_f = angle(exp(sqrt(-1).*S_tmp_pad)./exp(sqrt(-1).*S_tmp_pad_f));
		fprintf('high-pass\n');
	else
		fprintf('low-pass\n');
	end
	S_tmp_f = S_tmp_pad_f(N_pad(1)+1:end-N_pad(1), ...
                       	  N_pad(2)+1:end-N_pad(2), ...
                       	  N_pad(3)+1:end-N_pad(3));
	S_tmp_f(~SM_brain) = 0;
	S_phu_f(:, :, :, idx) = S_tmp_f;
end
%save_series([fname_in '_mag'], [fname_out '_phu'], single(S_phu), []);
%save_series([fname_in '_mag'], [fname_out '_phu_filt'], single(S_phu_f), []);
clear S_mag_in S_pha_in S_phu S_tmp_pad S_tmp_pad_f S_tmp_f;

% Reconstruct unwrapped phase by assuming a linear phase evolution: phi = k*TE + d
K_est = zeros(size(S_phu_f(:, :, :, 1)));
D_est = zeros(size(S_phu_f(:, :, :, 1)));
SoS_est = NaN(size(S_phu_f(:, :, :, 1)));
for z=1:size(K_est, 3)
    fprintf('%d ...', z);
    for y=1:size(K_est, 2)
        for x=1:size(K_est, 1)
			if SM_brain(x, y, z)
            	I = reshape(S_phu_f(x, y, z, :), N_te, 1); % was unwrap
            	if ~sum(isnan(I))
					P = quantreg(T(:), I(:), .5, 1);
					K_est(x, y, z) = P(1);
					D_est(x, y, z) = P(2);
					SoS_est(x, y, z) = sum(abs(I(:)-polyval(P, T(:))));
            	end
			end
        end
    end
end
fprintf('done.\n');

% Save result
S_f0 = K_est./(2*pi);
S_tmp = SoS_est(SM_brain);
S_tmp(isnan(S_tmp)) = [];
SM_975 = SoS_est < quantile(S_tmp(:), .975);
S_f0_975 = S_f0;
S_f0_975(~SM_975) = NaN;
SM_95 = SoS_est < quantile(S_tmp(:), .95);
S_f0_95 = S_f0;
S_f0_95(~SM_95) = NaN;
if expanded
	save_series([fname_in '_mag'], [fname_out '_f0'], ...
				single(cat(4, S_f0(:, :, 1:end-1), S_f0_975(:, :, 1:end-1), ...
							  S_f0_95(:, :, 1:end-1), SoS_est(:, :, 1:end-1))), []);
else
	save_series([fname_in '_mag'], [fname_out '_f0'], ...
				single(cat(4, S_f0, S_f0_975, S_f0_95, SoS_est)), []);
end

