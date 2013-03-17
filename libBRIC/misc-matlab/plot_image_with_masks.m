function [] = plot_image_with_masks(I, ti, M1, M2, M3, show_colorbar, I_max)
% Plots image I with masks M1, M2, M3 if given.
% INPUTS: I - 2D RGB or grayscale image
%         ti - image title
%         M1, M2, M3 - 2D binary masks
%         show_colorbar - flag indicating if colorbar should appear
%         I_max - maximum image intensity
%
if isempty(M1)
    nbit = 0;
    I_mask = zeros(size(I), 'uint8');
else
    nbit = 1; % bit number of respective mask
    I_mask = uint8(M1);
end
if ~isempty(M2)
    nbit = 2;
    I_mask = uint8(M2)*nbit + I_mask;
end
if ~isempty(M3)
    nbit = 4;
    I_mask = uint8(M3)*nbit + I_mask;
end
map = [ [0 0 0]; hsv(2^nbit-1) ]; % map contains all colors of rgb image
[I_rgb, tmp2] = level2rgb(I_mask, map); clear tmp2;
I_rgb = overlay_mask(I, I_rgb, I_max); % Assigns each level in I_ind a
                                       % different color
plot_image(I_rgb, ti, []);
if show_colorbar
	colormap(map); % Need to set colormap so that colorbar shows correct colours
	colorbar('SouthOutside');
end


