function [] = save_series(base_name, new_name, new_series, zidx)
% Load NIFTI or Analyze file (.nii.gz, .nii, or .hdr/.img)
% INPUTS: base_name - base image from which the header information (orientation)
%                     is extracted to create the header of the new file.
%         new_name - name of the NIFTI or Analyze file that should be created
%                    (If the file extension is not specified it is '.nii.gz')
%         new_series - the 3D, 4D, ... volume that should be saved in
%                      NIFTI or Analyze format. The output data type
%                      will be same as 'class(new_series)', except that logical
%                      masks will be saved as 'uint8'.
%         zidx - either [] if the output volume should have the same number of 
%                slices as size(new_series, 3) or the slice indices as
%                specified with 'load_series()'.
% EXAMPLES:
% addpath('NIFTI/');
% SM = load_series('1/13449/GRE', []) > 200; % binarize volume
% save_series('1/13449/GRE', '1/13449/GRE_bin', SM, []);
%
% SM = load_series('1/13449/GRE', 4:6) > 200; % binarize only slice 4:6
% save_series('1/13449/GRE', '1/13449/GRE_bin', SM, 4:6);
%

% Load the header of the base image
NII = load_series(base_name, 0);
if NII.hdr.dime.dim(1) ~= length(size(new_series))
    warning('save_series:dimerror', ...
            'Base volume and new volume have different dimensionality!');
end
% Add the image data
if isempty(zidx)
    if islogical(new_series)
        NII.img = uint8(new_series);
    else
        NII.img = new_series;
    end
else
    if islogical(new_series)
        NII.img = zeros(NII.hdr.dime.dim(2:NII.hdr.dime.dim(1)+1), 'uint8');
        NII.img(:, :, zidx, :, :) = uint8(new_series);
    else
        NII.img = zeros(NII.hdr.dime.dim(2:NII.hdr.dime.dim(1)+1), class(new_series));
        NII.img(:, :, zidx, :, :) = new_series;
    end
end

% Call make_nii to get correct datatype, bitpix value, ...
NII_tmp = make_nii(NII.img);
NII.hdr.dime.dim = NII_tmp.hdr.dime.dim;
NII.hdr.dime.datatype = NII_tmp.hdr.dime.datatype;
NII.hdr.dime.bitpix = NII_tmp.hdr.dime.bitpix;
NII.hdr.dime.glmax = NII_tmp.hdr.dime.glmax;
NII.hdr.dime.glmin = NII_tmp.hdr.dime.glmin;
NII.hdr.dime.cal_max = quantile(NII.img(:), .95);
NII.hdr.dime.cal_min = quantile(NII.img(:), .05);

% Add string about creator
NII.hdr.hist.descrip = 'Created by save_series()';

% Save series raw
save_series_raw(NII, new_name);
