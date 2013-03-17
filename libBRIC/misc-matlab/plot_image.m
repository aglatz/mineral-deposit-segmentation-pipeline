function [] = plot_image(I, ti, I_max)
% Plots the image I.
% INPUTS: I - 2D RGB or grayscale MRI image
%         ti - title
%         I_max - maximum image intensity
%
if size(I, 3) == 1
    if isempty(I_max)
        I_max = max(I(:));
    end
    [U, map] = gray2ind(mat2gray(I, double([0 I_max])), double(ceil(I_max)));
    U = ind2rgb(rot90(U), map);
else
    U = zeros(size(I, 2), size(I, 1), size(I, 3), class(I));
    for i = 1:size(I, 3)
        U(:,:,i) = rot90(I(:,:,i));
    end
end
image(U);
axis image;
axis off;
if ti; title(ti); end

