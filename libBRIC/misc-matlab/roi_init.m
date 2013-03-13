function [Roi] = roi_init(S)
% Takes a mask (logical up to uint16) of one continguos structure
% along the slice direction (e.g. basal ganglia) and calculates
% the slice indexes covered by the mask. Additonal processing is
% done to trim slices where the mask is false or 0 and to find
% the rectange encompassing the whole structure after maximum
% intensity projection along the slice direction and coversion to
% logical values.
% INPUTS: S - mask (logical up to uint16) of a continguous structure
%             along the slice direction.
% RETURNS: Roi - A structure containing the slice indices covered by
%                the mask, the trimmed mask, and the rectangle
%                that covers the maximum intensity projection of the
%                mask.
%
nslices = size(S, 3);
Roi.Num = zeros(nslices, 2);
Roi.Int = zeros(size(S), class(S));
mask_idx = 1;
for slice_idx = 1:nslices
    tmp = sum(S(:, :, slice_idx), 2);
    if sum(tmp) % Is there any basal ganglia structure in this slice?
        Roi.Num(mask_idx, 1) = slice_idx; % Nifti slice numbering
        Roi.Num(mask_idx, 2) = nslices-slice_idx+1; % Analyse slice numbering

		Roi.Int(:, :, mask_idx) = uint16(S(:, :, slice_idx));

        mask_idx = mask_idx + 1;
    end
end
Roi.Num = Roi.Num(1:mask_idx-1, :); % Trim away additional zeros
Roi.Int = Roi.Int(:, :, 1:mask_idx-1); % Trim away additional zeros

% Add information about the largest rectangle in X-Y-plane,
% which encompasses all BG structures in all slices.
[X, Y] = roi_maxrect(Roi, 0);
Roi.Maxrect.X = X;
Roi.Maxrect.Y = Y;
Roi.Maxrect.Mask = false(size(S(:, :, 1)));
Roi.Maxrect.Mask(X, Y) = ~Roi.Maxrect.Mask(X, Y);


%-------------------------------------------------------------------------
function [X_maxrect Y_maxrect] = roi_maxrect(Roi, perim_frac)
Size = roi_XYsize(Roi); % max X or Y extension
idx_x = [Size(1) 0];
idx_y = [Size(2) 0];
for sliceidx = 1:roi_nslices(Roi)
    [X Y] = roi_rect(Roi, sliceidx, perim_frac);
    idx_x(1) = min([ idx_x(1) X ]);
    idx_x(2) = max([ idx_x(2) X ]);
  
    idx_y(1) = min([ idx_y(1) Y ]);
    idx_y(2) = max([ idx_y(2) Y ]);
end
X_maxrect = idx_x(1):idx_x(2);
Y_maxrect = idx_y(1):idx_y(2);


%-------------------------------------------------------------------------
function [X, Y] = roi_rect(Roi, sliceidx, perim_frac)
Mask = roi_mask(Roi, sliceidx);

tmp = sum(Mask, 2);
X_range = [find(tmp>0, 1, 'first') find(tmp>0, 1, 'last')];

tmp = sum(Mask, 1);
Y_range = [find(tmp>0, 1, 'first') find(tmp>0, 1, 'last')];

% Add additional pixels to the ROI...
xperim = ceil((X_range(2) - X_range(1))/2 * perim_frac);
if X_range(1)-xperim > 0;               X_range(1) = X_range(1)-xperim; end;
if X_range(2)+xperim <= size(Mask, 1);	X_range(2) = X_range(2)+xperim; end;
yperim = ceil((Y_range(2) - Y_range(1))/2 * perim_frac);
if Y_range(1)-yperim > 0;               Y_range(1) = Y_range(1)-yperim; end;
if Y_range(2)+yperim <= size(Mask, 2);	Y_range(2) = Y_range(2)+yperim; end;

X = X_range(1):X_range(2);
Y = Y_range(1):Y_range(2);

% % Uncomment for testing
% figure;
% Mask1 = false(size(Mask));
% Mask1(X,Y) = ~Mask1(X,Y);
% Mask1 = bwperim(Mask1);
% rgb = imoverlay(mat2gray(Mask), Mask1, [1 0 0]);
% imagesc(rgb);
% input('enter');


