addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/LIBRA');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

S = load_series('${ARG_0}', []);
% Sum Z
S_z = sum(S, 3);
% Sum X
S_zx = sum(S_z, 1);
t = quantile(S_zx, ${ARG_1});
Y = find(S_zx > 0 & S_zx < t);
% Sum Y
S_zy = sum(S_z, 2); 
t = quantile(S_zy, ${ARG_1});
X = find(S_zy > 0 & S_zy < t);

% Create mask
SM = false(size(S));
SM(X, Y, :) = true;
save_series('${ARG_0}', '${ARG_2}', SM, []);

N = mscalelogist(double(S(SM)))/.655;

% Save result
fd = fopen(['${ARG_2}' '.txt'], 'w');
fprintf(fd, '%d', round(N));
fclose(fd);

fprintf('%0.3f\n', N);
