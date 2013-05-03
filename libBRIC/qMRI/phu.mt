addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

% Input data
fname_in = '${ARG_0}';
N_te = ${ARG_1};
te = ${ARG_2}; % in us
te2 = ${ARG_3}; % in us
fname_out = '${ARG_4}';

% Calculate TEs
if N_te > 1
    T = (te:(te2-te):te+(te2-te)*(N_te-1))';
else
    T = te;
end

% Load input data
S_abs = load_series([fname_in '_mag'], []);
S_pha = load_series([fname_in '_pha'], []);


% Reconstruct unwrapped phase by assuming a linear phase evolution: phi = k*TE + d
% Calculate first estimate for: k=deta_phi/delta_TE; d=phi-k*TE
N_diff = N_te-1;
K = zeros(size(S_abs(:, :, :, 1), 1), size(S_abs(:, :, :, 1), 2), ...
          size(S_abs(:, :, :, 1), 3), N_diff);
for idx = 1:N_diff
    S_diff = ((S_abs(:, :, :, idx+1) .* exp(sqrt(-1)*S_pha(:, :, :, idx+1))) ./ ...
              (S_abs(:, :, :, idx) .* exp(sqrt(-1)*S_pha(:, :, :, idx))));
    
    % remove NaNs and Infitity
    SM = S_abs(:, :, :, idx) == 0;
    S_diff(SM) = zeros(sum(SM(:)), 1);
    
	K(:, :, :, idx) = angle(S_diff) ./ (T(2)-T(1)); % rad/us
end
K_avg = mean(K, 4); % in rad/us
clear K;

% Field map: dphi=-gyro*dTE*B0
F0 = (K_avg ./ (2*pi)) .* 1E6; % Hz
save_series([fname_in '_mag'], [fname_out '_F0'], single(F0), []);

% Phase offset
D = zeros(size(S_abs(:, :, :, 1), 1), size(S_abs(:, :, :, 1), 2), ...
          size(S_abs(:, :, :, 1), 3), N_te);
for idx = 1:N_te
    S_diff = ((S_abs(:, :, :, idx) .* exp(sqrt(-1) .* S_pha(:, :, :, idx))) ./ ...
              (S_abs(:, :, :, idx) .* exp(sqrt(-1) .* K_avg .* T(idx) .* ones(size(S_pha(:, :, :, idx))))));
	D(:, :, :, idx) = angle(S_diff);
end
D_avg = mean(D, 4); % in rad
clear D;

% Calculate phase jumps
S_diff = zeros(size(S_pha), class(S_pha));
for idx = 1:N_te
    S_tmp = single(D_avg + K_avg * T(idx));
    S_diff(:, :, :, idx) = S_tmp - S_pha(:, :, :, idx);
end
clear K_avg S_tmp D_avg;

% Unwrap field map
S_phu = S_pha + round(S_diff./(2*pi)) * (2*pi);

% Fit linear phase evolution and denoisify field map
K_est = zeros(size(S_phu(:, :, :, 1)), 'single');
D_est = zeros(size(S_phu(:, :, :, 1)), 'single');
for z=1:size(K_est, 3)
    fprintf('%d ...', z);
    for y=1:size(K_est, 2)
        for x=1:size(K_est, 1)
            I = reshape(S_phu(x, y, z, :), N_te, 1);
            if ~sum(isnan(I))
                P = polyfit(T, I, 1);
                K_est(x, y, z) = P(1);
                D_est(x, y, z) = P(2);
            end
        end
    end
end
fprintf('done.\n');

% Save result
for idx = 1:N_te
    S_phu(:, :, :, idx) = single(D_est + K_est * T(idx));
end
save_series([fname_in '_mag'], [fname_out '_phu'], S_phu, []);

