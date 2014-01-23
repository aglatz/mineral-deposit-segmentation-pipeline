% Internal function

function [I_level_rgb, map] = level2rgb(I_level, map)
% Sanity check
if max(I_level(:)) > intmax('uint8') && min(I_level(:)) < intmin('uint8')
    warning('Down-casting I_level to uint8 will change original data!');
end

% I_ind HAS to be uint8 otherwise ind2rgb() shows a different behaviour
if min(I_level(:)) ~= 0
    I_tmp = uint8(I_level);
    I_ind = uint8(I_tmp - min(I_tmp(:))*ones(size(I_level), class(I_tmp)));
else
    I_ind = uint8(I_level);
end

% We need size(map, 1) == max(I_ind(:))+1
if isempty(map)
    col_levels = double(max(I_ind(:)));
    map = [ [0 0 0]; hsv(col_levels) ]; % Add background
else
    if size(map, 1) < max(I_ind(:))+1
        warning('RGB map too short. Different indexed pixels will have the same color.');
    end
end

% Convert
I_level_rgb = ind2rgb(I_ind, map);

