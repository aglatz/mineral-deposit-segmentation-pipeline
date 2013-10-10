addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

% Read volume
S = load_series('${ARG_0}', []);

% Read z coordinates
fd = fopen('${ARG_1}', 'r');
Idx = fscanf(fd, '%d %d');
fclose(fd);

% Save trucated volume
save_series('${ARG_0}', '${ARG_2}', S(:, :, Idx(1):Idx(2)), []);
