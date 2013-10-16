function [S_r2smap, S_r2smap_sd, S_s0map, S_s0map_sd, S_csqmap, S_r2slog] = ...
                recon_r2smap_lmf(S_abs, N, T, threshold, verbose)
% This function fits mono-exponential functions to the multi-echo
% input data and estimates the transversal realaxation time for every voxel.
% Inputs: S_abs - 4D volume where the 4th dimension is the time
%         N - Number of echos to use for the fit
%         T - Echo times in us
%         threshold - Skips voxels with first echo data below this threshold
%         verbose - flag indicating if function should be more verbose
% Outputs: S_r2smap - Estimated transversal relaxation times R2
%          S_r2smap_sd - Estimated standard error of R2
%          S_s0map - Estimated signal S(t=0s)
%          S_s0map_sd - Estimated standard error of S(t=0s)
%          

S_r2smap = zeros(size(S_abs(:, :, :, 1)), 'single');
S_r2smap_sd = zeros(size(S_abs(:, :, :, 1)), 'single');
S_s0map = zeros(size(S_abs(:, :, :, 1)), 'single');
S_s0map_sd = zeros(size(S_abs(:, :, :, 1)), 'single');
S_csqmap = zeros(size(S_abs(:, :, :, 1)), 'single');
S_r2slog = zeros(size(S_abs(:, :, :, 1)), 'uint8');

if ~isempty(verbose) && verbose
    fprintf('Reconstructing r2s map...');
end

T = double(T);
N = double(N);

for s = 1:size(S_abs, 3)
    if ~isempty(verbose) && verbose
        fprintf('%d...', s);
    end
	
    % Per-voxel fitting
    for j = 1:size(S_abs, 1)
        for m = 1:size(S_abs, 2)
            % Per-voxel time signal and noise
            I = double(reshape(S_abs(j, m, s, :), size(S_abs, 4), 1));
            
            % 3SD above noise is good to go...
            M = I > 3*N;
            if I(1) > threshold && sum(M) > length(I)*0.5; 
                % Estimate starting values for the fitting procedure
                P_linfit = polyfit(T(M), log(I(M)), 1);
                S_r2slog(j, m, s) = 1;
                try
                    fun = @(c) mrqcof(T(M), I(M), I(M)./N(M), c, 0);
                    [res, ssq, cnt] = LMFnlsq(fun, [exp(P_linfit(2)) -P_linfit(1)]); %, 'Display', 1);
                    S_s0map(j, m, s) = res(1);
                    S_r2smap(j, m, s) = res(2);
                catch
                    S_r2slog(j, m, s) = 2;
                end
                if S_r2slog(j, m, s) == 1
                    % Ok if there is something inside...
                    if S_s0map(j, m, s) > 0 && S_r2smap(j, m, s) > 0 && cnt > 0
                        % Sum of squares
                        [dy, alpha] = mrqcof(T(M), I(M), I(M)./N(M), ...
                            [S_s0map(j, m, s) S_r2smap(j, m, s)], 0);
                        S_csqmap(j, m, s) = sum(dy.^2);

                        % Estimated error of parameters
                        ma = 2;
                        covar = zeros(ma);
                        unitmat = eye(ma);
                        for i=1:ma
                            inhom = unitmat(:, i);
                            augmented = [alpha inhom];
                            cr = rref(augmented);
                            covar(:, i) = cr(:, ma+1);
                        end
                        S_s0map_sd(j, m, s) = sqrt(covar(1, 1));
                        S_r2smap_sd(j, m, s) = sqrt(covar(2, 2));

                        % Success
                        S_r2slog(j, m, s) = 4;
                    else
                        S_r2slog(j, m, s) = 3;
                    end
                end
            end
        end
    end
end

if ~isempty(verbose) && verbose
    fprintf('done.\n');
end


%--------------------------------------------------------------------------
function [dy, alpha] = mrqcof(x, y, sig, func_a, func_c)
ndata = numel(x);
ma = numel(func_a);
alpha = zeros(ma, ma);

[ymod, dyda] = func_monoexp_const(x, func_a, func_c);

dy = (ymod - y) .* (sig./max(sig));
for k = 1:ndata
    for i = 1:ma
        for j = 1:ma
            alpha(i,j) = alpha(i,j) + (1./(sig(k).^2))*dyda(i,k)*dyda(j,k);
        end
    end
end


%--------------------------------------------------------------------------
function [ymod, dyda] = func_monoexp_const(x, a, c)
t = x;
ma = numel(a);
ndata = numel(x);

dyda = zeros(ma, ndata);
ymod = zeros(ndata, 1);
for j=1:ndata
    dyda(1,j) = exp(-a(2) * t(j));
    dyda(2,j) = a(1) * exp(-a(2) * t(j)) .* -t(j);
    ymod(j) = a(1) * exp(-a(2) * t(j)) + c;
end
