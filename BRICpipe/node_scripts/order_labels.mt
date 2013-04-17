addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

S = load_series('${ARG_0}', []);
S_lab = load_series('${ARG_1}', []);
Lab = unique(S_lab(:))';
Lab = Lab(Lab > 0 & Lab < 256);
N_lab = length(Lab);
Means = zeros(N_lab, 1);
for idx = 1:N_lab
    S_tmp = S .* cast(S_lab == Lab(idx), class(S));
    S_tmp(~S_tmp) = [];
    Means(idx) = mean(S_tmp);
end
[tmp, idx_sorted] = sort(Means, 1);
fd = fopen('${ARG_2}', 'w');
for idx = 1:N_lab
    fprintf(fd, '%d ', Lab(idx_sorted(idx)));
end
fclose(fd);

