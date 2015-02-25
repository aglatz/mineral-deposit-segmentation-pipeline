addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

S = load_series('${ARG_0}', []);
dz_mm = 0;
switch size(S, 4)
    case 1, S = S;
    case 2, S = S(:, :, :, 2); % t2w
    otherwise, S = S(:, :, :, 3); % 3th echotime - less sus. artefacts
end
N_slice = size(S, 3);
I_med = NaN(N_slice, 1);
for idx_slice = 1:N_slice
    S_tmp = S(:, :, idx_slice);
    I_tmp = S_tmp(:);
    I_med(idx_slice) = median(double(I_tmp(I_tmp > 0)));
end
Idx = find(I_med > quantile(I_med(I_med > 0), 0));
if dz_mm > 0
    NII = load_series('${ARG_0}', 0);
    dz = round(dz_mm/NII.hdr.dime.pixdim(4));
else
    dz = 0;
end
idx_min = min(Idx)+dz;
idx_max = max(Idx);

% Save for later
fd = fopen('${ARG_1}.txt', 'w');
fprintf(fd, '%d %d', idx_min, idx_max);
fclose(fd);

% Save trucated volume
save_series('${ARG_0}', '${ARG_1}', S(:, :, idx_min:idx_max), []);

