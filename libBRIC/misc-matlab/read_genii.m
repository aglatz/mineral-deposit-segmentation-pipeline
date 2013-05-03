function [raw_abs, raw_pha] = read_genii(fname, T, use_reaima, xyres)
% Reads multi-echo MRI data from a GE scanner and converts it
% into 4D NIFTI magnitude and phase volumes, where the phase
% ranges from 0 to 2*pi). The input volumes' file names should be
% in the format <path+basename>_<echotime in us>.<extension>.
% Inputs: fname - path+basename of input volumes
%         T - echo times of input volumes
%         use_reaima - a flag indicating if the final magnitude and
%                      phase volumes should be calculated from the
%                      real and imaginary input volumes
%         xyres - xy resolution of final 4D volumes or empty
%
Hdr = load_series(sprintf('%s_%d', fname, T(1)), 0);
if ~isempty(xyres)
    sz = [xyres Hdr.hdr.dime.dim(4) numel(T)];
else
	sz = [Hdr.hdr.dime.dim(2:4) numel(T)];
end
raw_abs = zeros(sz, 'single');
raw_pha = zeros(sz, 'single');

for idx = 1:numel(T)
    I_tmp = single(load_series(sprintf('%s_%d', fname, T(idx)), []));
    if use_reaima
        I_cpx = complex(I_tmp(:, :, :, 3), I_tmp(:, :, :, 4));
        for slice = 1:size(I_cpx, 3)
            if ~mod(slice, 2)
                I_cpx(:, :, slice) = -1 * I_cpx(:, :, slice); % prelude: 0..2*pi
            end
        end
    else
        if ndims(I_tmp) > 3
            I_pha = (I_tmp(:, :, :, 2)+pi*1000)/1000; % prelude: 0..2*pi
        else
            I_pha = zeros(size(I_tmp), 'single');
        end
        I_cpx = I_tmp(:, :, :, 1) .* exp(sqrt(-1) .* I_pha);
    end
    if ~isempty(xyres)
        I_cpx = fftInterpolate(I_cpx, [xyres size(I_cpx, 3)]);
    end
    raw_abs(:, :, :, idx) = abs(I_cpx);
    raw_pha(:, :, :, idx) = angle(I_cpx);
end
