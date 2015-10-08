function [gehdr, out_abs, out_pha] = read_pfile(pfile, recon)
% Reads GE pfile data (kspace) and transforms them into image space.
% INPUT: pfile - pfile name
%        recon - if 'true' it combines the images from a coil array (see
%                below)
% OUTPUT: gehdr - pfile header as a struct
%         out_abs - reconstructed magnitude images
%         out_pha - reconstructed phase images
%

% TODO make them a varin
x_idx = [];
y_idx = [];
slices = [];
echos = [];
SM_noise = [];

fd = fopen(pfile, 'r', 'ieee-le');

% Read header
gehdr = read_gehdr(fd);
hdr = get_pfile_hdr(gehdr);

if isempty(x_idx)
    x_size = hdr.nframe-1;
    x_idx = 1:x_size;
else
    x_size = length(x_idx);
end

if isempty(y_idx)
    y_size = hdr.nsamp;
    y_idx = 1:y_size;
else
    y_size = length(y_idx);
end

if isempty(slices)
    z_idx = 1:hdr.nslice;
    slices = z_idx;
    z_size = hdr.nslice;
else
    z_size = length(slices);
    z_idx = 1:z_size;
end

if isempty(echos)
    echos_idx = 1:hdr.necho;
    echos = echos_idx;
    necho = hdr.necho;
else
    necho = length(echos);
    echos_idx = 1:necho;
end

% Read data
raw_data = zeros(x_size, y_size, z_size, necho, hdr.nrecv);
for slice_idx = z_idx
    for echo_idx = echos_idx
        slice_data = read_pfile_data(fd, hdr, echos(echo_idx), slices(slice_idx));
        raw_data(:, :, slice_idx, echo_idx, :) = slice_data(x_idx, y_idx, :);
    end
end

fclose(fd);

weight = ones(hdr.nrecv, 1);
if ~isempty(SM_noise)
    N_mean = NaN(hdr.nrecv, 1);
    N_std = NaN(hdr.nrecv, 1);
    for recv_idx = 1:hdr.nrecv
        S_tmp = abs(raw_data(:, :, :, 1, recv_idx));
        N_mean(recv_idx) = mloclogist(S_tmp(SM_noise));
        N_std(recv_idx) = mscalelogist(S_tmp(SM_noise));
        fprintf('recv:%d N_mean=%0.2f N_std=%0.2f\n', ...
                recv_idx, N_mean(recv_idx), N_std(recv_idx));
    end
    weight = (sum(N_std)/hdr.nrecv)./N_std
end

if recon
    % See: Neuroimage, 2010 Jul 1; 51(3):1089-97
    out_data = zeros(x_size, y_size, z_size, necho-1);
    for echo_idx = 2:necho
        S_tmp = zeros(size(raw_data(:, :, :, 1, 1)));
        for recv_idx = 1:hdr.nrecv;
            S_tmp = S_tmp + raw_data(:, :, :, echo_idx, recv_idx) .* ...
                            conj(raw_data(:, :, :, 1, recv_idx)) .* ...
                            weight(recv_idx);
        end
        out_data(:, :, :, echo_idx-1) = S_tmp ./ hdr.nrecv;
    end
    out_abs = abs(out_data);
    out_pha = angle(out_data);
else
    out_abs = abs(raw_data);
    out_pha = angle(raw_data);
end


% function [out_abs, out_pha] = do_recon(in, recon)
% if recon
%     in_abs = abs(in);
%     out_pha = angle(in);
% 	idx_x = get_center_idx(size(out_pha, 1));
%     idx_y = get_center_idx(size(out_pha, 2));
%     for rcv = 1:size(out_pha, 3)
%         px_center = out_pha(idx_x, idx_y, rcv);
%         out_pha(:, :, rcv) = angle( exp(sqrt(-1) * out_pha(:, :, rcv)) ./ ...
%                              exp(sqrt(-1) * mean(px_center(:)) * ...
%                              ones(size(out_pha(:, :, rcv)), class(out_pha))) );
%     end
%     in_corr = in_abs .* exp(sqrt(-1) .* out_pha);
% 	out_abs = sqrt(sum(in_abs.^2, 3));
%     out_abs = out_abs ./ max(out_abs(:)) * 4095; % Scale as GE possibly does it...
%     out_pha = angle(sum(in_corr, 3) ./ sum(in_abs, 3)) + pi; % range 0 to 2*pi for fsl prelude
% else
%     out_abs = abs(in);
%     out_pha = angle(in) + pi; % range 0 to 2*pi for fsl prelude
% end
% 
% 
% function [idx] = get_center_idx(idx_max)
% if mod(idx_max, 2) == 0
% 	idx = [idx_max/2, idx_max/2 + 1];
% else
% 	idx = ceil(idx_max/2);
% end

