clear all; close all;

% addpath ~/tmp/mineral/NIFTI
% addpath ~/tmp/mineral/LMFnlsq
% addpath ~/tmp/mineral-deposit-segmentation-pipeline/libBRIC/misc-matlab

addpath /ISIS/proc1/aglatz/mineral-deposit-segmentation-pipeline/libBRIC/qMRI/LMFnlsq/
addpath /ISIS/proc1/aglatz/mineral/NIFTI/
addpath /ISIS/proc1/aglatz/mineral-deposit-segmentation-pipeline/libBRIC/misc-matlab/

cdir = '.';
S = double(load_series(fullfile(cdir, 'S_reg'), []));
G = double(load_series(fullfile(cdir, 'S_phu_G'), []));
for i=1:8
    S_tmp = S(:, :, i);
    S_tmp(isnan(G)) = 0;
    S(:, :, i) = S_tmp;
end

T = [4 7 11.3 12 17 24.7 35 42]*1e-3
T = T(:);

[S_r2smap, S_r2smap_sd, S_s0map, S_s0map_sd, S_csqmap, S_r2slog] = ...
	recon_r2smap_lmf(S, ones(size(S, 4), 1)*10, T, 10, G);

save_series(fullfile(cdir, 'S'), fullfile(cdir, 'R2s'), S_r2smap, []);
save_series(fullfile(cdir, 'S'), fullfile(cdir, 'S0'), S_s0map, []);
save_series(fullfile(cdir, 'S'), fullfile(cdir, 'R2s_csq'), S_csqmap, []);
save_series(fullfile(cdir, 'S'), fullfile(cdir, 'R2s_log'), S_r2slog, []);
