addpath('${SCRIPT_DIR}/NIFTI');
addpath('${SCRIPT_DIR}/libBRIC/misc-matlab');

SM_basgang = load_series('${ARG_0}', []);
S_roi = zeros(size(SM_basgang), 'uint8');
for i = ...
[13 52 11 12 50 51 10 49 16] % FSL FIRST labels: regions to dilate
							 % come first!
							 % [13 12 11 10] % L
							 % [52 51 50 49] % R
	if i == 13 || i == 52
		SE = strel('disk', 6);
		SM_tmp = logical(imdilate(SM_basgang == i, SE));
		% New regions are called 14 and 55
		if i == 13
			i_new = 14;
		else
			i_new = 55;
		end
		S_roi(SM_tmp) = i_new;
	end
	SM_tmp = SM_basgang == i;
	S_roi(SM_tmp) = i;
end
% Remove thalamus and brainstem
for i = ...
[10 49 16]
	SM_tmp = S_roi == i;
	S_roi(SM_tmp) = 0;
end
save_series('${ARG_0}', '${ARG_1}', S_roi, []);

