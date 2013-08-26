function [S] = load_series(basename, zidx)
% Load NIFTI or Analyze file (.nii.gz, .nii, or .hdr/.img)
% INPUTS: basename - input file name
%         zidx - (i) a vector containing valid slice numbers, or
%                (ii) 0, which means that just the header should be returned, or
%                (iii) [], which means that all slices should be returned
% RETURNS: S - only the image data without the header, or
%              only the header without image data
% EXAMPLE:
% addpath('NIFTI/');
% S = load_series('1/13449/GRE', []);
%
[path, name, ext] = fileparts(basename);
if strcmp(ext, '.gz')
    [tmp1, name, tmp2] = fileparts(name);
end
if ~isempty(path)
    basename = fullfile(path, name);
else
    basename = name;
end
name = [basename '.hdr'];
fd = fopen(name, 'r');
if fd < 0
    name = [basename '.nii.gz'];
    fd = fopen(name, 'r');
    if fd < 0
        name = [basename '.nii'];
        fd = fopen(name, 'r');
        if fd < 0
            error(['Unrecognized file type: ' basename]);
        else
            fclose(fd);
            [S] = load_series_core(name, zidx);
        end
    else
        fclose(fd);
        [S] = load_nifti_gz_series(name, zidx);
    end
else
    fclose(fd);
    [S] = load_series_core(name, zidx);
end


function [S] = load_nifti_gz_series(name, zidx)
% We prefer files in NIFTI_GZ format since this saves a lot of
% disk space and we have a 'unzip' command. The MATLAB gunzip 
% command is provided by a Java class, so this command requires 'jvm'.
tmp_dir = tempname('.');
[ret, msg] = mkdir(tmp_dir);
if ~strcmp(msg, '')
    error(['Problems creating temporary directory: ' tmp_dir]);
end
try
    % Extract ...
    file_path = char(gunzip(name, tmp_dir));
    % Load ...
    [S] = load_series_core(file_path, zidx);
    % Cleanup...
    rmdir(tmp_dir, 's');
catch
    error(['Problems reading file: ' name]);
end


function [S] = load_series_core(name, zidx)
%
if length(zidx) == 1 && ~zidx
    NII = load_untouch_nii(name);
    S = rmfield(NII, 'img'); % we just want the header
else
    NII = load_untouch_nii(name);
    if NII.hdr.dime.dim(1) > 5 % 5D is maximum
        error('load_series:load_series_core:maxdim', ...
              'Input volume has more than 5 dimensions!');
    end
    S = NII.img;
    if ~isempty(zidx)
        S = S(:, :, zidx, :, :);
    end
end



