% Internal function

function [I_gray_rgb] = overlay_mask(I_gray, I_seg_rgb, I_max)
% Generate binary mask from rgb mask
M_seg = sum(I_seg_rgb, 3) ~= 0;

[I_tmp, Map_gray] = gray2ind(mat2gray(I_gray, double([0 I_max])), double(ceil(I_max)));

% Convert into rgb image so we are independent of colormap when
% plotting the image and add mask colours
I_gray_rgb = ind2rgb(I_tmp, Map_gray);
for i=1:size(I_gray_rgb, 3);
    I_tmp_gray = I_gray_rgb(:, :, i);
    I_tmp_seg = I_seg_rgb(:, :, i);
    I_tmp_gray(M_seg) = I_tmp_seg(M_seg);
    I_gray_rgb(:, :, i) = I_tmp_gray;
end

