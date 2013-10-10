addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

% Input data
fname_in = '${ARG_0}';
N_te = ${ARG_1};
te = ${ARG_2}; % in us
te2 = ${ARG_3}; % in us
bmask = '${ARG_5}';
fname_out = '${ARG_6}';

% Calculate TEs
if N_te > 1
    T = (te:(te2-te):te+(te2-te)*(N_te-1))';
else
    T = te;
end

% Load input data
S_abs = single(load_series([fname_in '_mag'], []));
S_pha = single(load_series([fname_in '_pha'], []));
SM_brain = logical(load_series(bmask, []));

% Remove low frequency variations
NII = load_series([fname_in '_mag'], 0);
F = NII.hdr.dime.pixdim(2:4);
for idx=1:N_te
	S_tmp = S_abs(:, :, :, idx) .* exp(sqrt(-1)*S_pha(:, :, :, idx));
	S_tmp_f = S_tmp ./ gauss3filter(S_tmp, ${ARG_4}, F);
	S_abs(:, :, :, idx) = abs(S_tmp_f) .* cast(SM_brain, class(S_tmp_f));
	S_pha(:, :, :, idx) = angle(S_tmp_f) .* cast(SM_brain, class(S_tmp_f));
end
save_series([fname_in '_mag'], [fname_out '_filt'], single(S_pha), []);

% Reconstruct unwrapped phase by assuming a linear phase evolution: phi = k*TE + d
K_est = zeros(size(S_pha(:, :, :, 1)), 'single');
D_est = zeros(size(S_pha(:, :, :, 1)), 'single');
SoS_est = zeros(size(S_pha(:, :, :, 1)), 'single');
for z=1:size(K_est, 3)
    fprintf('%d ...', z);
    for y=1:size(K_est, 2)
        for x=1:size(K_est, 1)
			if SM_brain(x, y, z)	
            	I = unwrap(reshape(S_pha(x, y, z, :), N_te, 1));
            	if ~sum(isnan(I))
                	P = polyfit(T, I, 1);
                	K_est(x, y, z) = P(1);
                	D_est(x, y, z) = P(2);
					SoS_est(x, y, z) = sum((I - polyval(P, T)).^2);
            	end
			end
        end
    end
end
fprintf('done.\n');

% Save result
S_phu = zeros(size(S_pha));
for idx = 1:N_te
    S_phu(:, :, :, idx) = single(D_est + K_est * T(idx));
end
save_series([fname_in '_mag'], [fname_out '_phu'], S_phu, []);
save_series([fname_in '_mag'], [fname_out '_sos'], SoS_est, []);
SM = SoS_est > quantile(SoS_est(SM_brain), .05) & SoS_est < quantile(SoS_est(SM_brain), .95);
save_series([fname_in '_mag'], [fname_out '_sos_bin'], SM, []);

