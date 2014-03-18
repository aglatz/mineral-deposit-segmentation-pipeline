addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

S_gre = double(load_series('${ARG_0}', [])); % Magnitude
S_gre = S_gre(:, :, :, ${ARG_1});
TE = ${ARG_2};
S_roi = load_series('${ARG_3}', []); % ROIs
Tmp = textscan('${ARG_4}', '%d,');
Idx_roi = Tmp{1}; Idx_roi = Idx_roi(:)';

% Single exponential fit
fprintf('Roi idxs:\n');
unique(S_roi(:))'
SM = false(size(S_roi));
for idx_roi = Idx_roi
	SM = SM | S_roi == idx_roi;
end
I_med = median(S_gre(SM));
fprintf('Roi:'); Idx_roi
fprintf('I_med=%0.2f\n', I_med);
R2s = -1/TE * log(S_gre./median(S_gre(SM)));

% Save results
[path, name, ext] = fileparts('${ARG_5}');
if isempty(path)
	path = '.';
end
if strcmp(ext, '.gz')
	[tmp, name, ext] = fileparts(name);
	ext = [ext '.gz'];
end
save_series('${ARG_0}', fullfile(path, name), single(R2s), []);

