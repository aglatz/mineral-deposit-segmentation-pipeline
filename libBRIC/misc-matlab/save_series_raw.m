function [] = save_series_raw(NII, basename)
% Save as NIFTI or Analyze file (.nii.gz, .nii, or .hdr/.img)
%
% This is an internal function and it's recommended to use
% the function 'save_series()' instead.
%
do_gzip = false;
[path, name, ext] = fileparts(basename);	
if ~isempty(ext)
    if strcmp(ext, '.gz')
        do_gzip = true;
        if ~isempty(path)
            name_ext = fullfile(path, name);
        else
            name_ext = name;
        end
    else
        do_gzip = false;
        name_ext = basename;
    end
else
    do_gzip = true;
    name_ext = [basename '.nii'];
end
if ~isfield(NII,'untouch') || NII.untouch == 0
	save_nii(NII, name_ext);
else
    save_untouch_nii(NII, name_ext);
end

try
    if do_gzip
        % Compress and cleanup - needs 'jvm'!
        gzip(name_ext);
        delete(name_ext);
    end
catch
    error(['Could not write file: ' name]);
end
