function [data] = read_pfile_data(fd, hdr, echo, slice)
% This script is based on the one from F. Frigo; see here:
% www.eng.mu.edu/frigof

% data_offset, point_size, nsamp, nframe, nrecv, necho, nslice
if echo > hdr.necho || slice > hdr.nslice
    error('Echo or slice number out of range');
end

% Init 3D matrix
data = zeros(hdr.nframe-1, hdr.nsamp, hdr.nrecv);

% Compute size (in bytes) of each frame, echo and slice.
% Here '2' counts real and imaginary value.
data_elements = (hdr.nsamp * 2) * (hdr.nframe - 1);	% number of data elements
                                                    % of a slice (no baseline view)
frame_size = (hdr.nsamp * 2) * hdr.point_size;      % size in bytes of a data frame
echo_size = frame_size * hdr.nframe;                % size in bytes of a slice
slice_size = echo_size * hdr.necho;                 % size in bytes of all echos of a slice
volume_size = slice_size * hdr.nslice;              % size in bytes of all echos and slices
                                                    % of a volume

for recv = 1 : hdr.nrecv
    % Compute offset in bytes to start of frame.  (skip baseline view)
    offset = hdr.data_offset ...
                + ((recv - 1) * volume_size) + ...
                + ((slice - 1) * slice_size) + ...
                + ((echo - 1) * echo_size) + ...
                + (frame_size);

    fseek(fd, offset, 'bof');

    % read data: element_size = 2 means 16 bit data,
    % element_size = 4 means EDR
    if hdr.point_size == 2
        raw_data = fread(fd, data_elements, 'integer*2');
    else
        raw_data = fread(fd, data_elements, 'integer*4');
    end

    for frame = 1 : (hdr.nframe-1)
        row_offset = (frame-1) * (hdr.nsamp*2);
        for samp = 1 : hdr.nsamp
            data(frame, samp, recv) =   raw_data( ((2*samp)-1) + row_offset ) + ...
                                        sqrt(-1)*raw_data( (2*samp) + row_offset );
        end                
    end
    
    % Each slice contains the raw data... (k-space!). The ifft2 along
    % z was calculated by the scanner.
    % Now we just need ifft2(). Use fftshift/ifftshift in to
    % handle the fact that kx,ky == 0 is in the middle of the matrix.
    % See e.g. blog.mshalin.com/2009/12/interplay-of-fft-ifft-fftshift-and_750.html
    % TODO solve the scaling issue?
    data(:, end:-1:1, recv) = fftshift( ifft2( ifftshift( data(:, :, recv) ) ) ) * 1024;
end