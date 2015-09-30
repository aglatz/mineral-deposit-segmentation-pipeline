clear all; close all;

% addpath ~/tmp/mineral/NIFTI
% addpath ~/tmp/mineral/LMFnlsq
% addpath ~/tmp/mineral-deposit-segmentation-pipeline/libBRIC/misc-matlab

addpath /ISIS/proc1/aglatz/mineral-deposit-segmentation-pipeline/libBRIC/qMRI/LMFnlsq/
addpath /ISIS/proc1/aglatz/mineral/NIFTI/
addpath /ISIS/proc1/aglatz/mineral-deposit-segmentation-pipeline/libBRIC/misc-matlab/

cdir = '.';
S = double(load_series(fullfile(cdir, 'S'), []));
T = [20 40 60 80]*1e-3
T = T(:);

[S_r2smap, S_r2smap_sd, S_s0map, S_s0map_sd, S_csqmap, S_r2slog] = ...
	recon_r2smap_lmf(S, ones(size(S, 4), 1), T, 10, 0);

save_series(fullfile(cdir, 'S'), fullfile(cdir, 'R2'), S_r2smap, []);
save_series(fullfile(cdir, 'S'), fullfile(cdir, 'S0'), S_s0map, []);
save_series(fullfile(cdir, 'S'), fullfile(cdir, 'R2_csq'), S_csqmap, []);
save_series(fullfile(cdir, 'S'), fullfile(cdir, 'R2_log'), S_r2slog, []);
