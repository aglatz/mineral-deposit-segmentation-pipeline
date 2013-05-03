function [] = save_series_kspace(S_abs, S_pha, fname, T)
% Save k-space in complex analyze format so it can be read by GUIBOLD.
% Inputs: S_abs - 4D magnitude data
%         S_pha - 4D phase data
%         fname - path and basename
%         T - echo times
%
for tpoint = 1:size(S_abs, 4)
    save_volume_kspace(S_abs(:, :, :, tpoint), S_pha(:, :, :, tpoint), ...
                       [fname '_' num2str(T(tpoint))]);
end


%--------------------------------------------------------------------------
function [] = save_volume_kspace(S_abs, S_pha, fname)
X = S_abs .* exp(sqrt(-1)*S_pha);

% Use fftshift/ifftshift in case we have to handle odd numbered frames.
X = fftshift( fftn( ifftshift( X ) ) );

XX = zeros(numel(X)*2, 1, 'single'); % "*2" bc 1pix = [real, imag]
XX(1:2:end) = reshape( single(real(X)), numel(X), 1 );
XX(2:2:end) = reshape( single(imag(X)), numel(X), 1 );

fd = fopen([fname '.img'], 'w', 'ieee-le');
if ~fd
    error(['Cannot open file: ' fname '.img']);
end
n = fwrite(fd, XX', 'single');
if n ~= length(XX)
    error(['Problems writing data to file: ' fname '.img']);
end
fclose(fd);
