addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/LIBRA');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

S = load_series('${ARG_0}', []);
try
	SM = logical(load_series('${ARG_1}', []));
catch
	f = str2double('${ARG_1}');
	S = S(:, :, :, 1);
	% Sum Z
	S_z = sum(S, 3);
	% Sum X
	S_zx = sum(S_z, 1);
	t = quantile(S_zx, f);
	Y = find(S_zx > 0 & S_zx < t);
	% Sum Y
	S_zy = sum(S_z, 2); 
	t = quantile(S_zy, f);
	X = find(S_zy > 0 & S_zy < t);

	% Create mask
	SM = false(size(S));
	SM(X, Y, :) = true;
	save_series('${ARG_0}', '${ARG_2}', SM, []);
end

N_mean = mloclogist(double(S(SM)))
N_std = mscalelogist(double(S(SM)))
N = N_std/.655

% Save result
fd = fopen(['${ARG_2}' '.txt'], 'w');
fprintf(fd, '%d', round(N));
fclose(fd);

