function [CC] = conncomp_list(L, L_iso, S_roi, S_weight, F, SM_hypo, SM_hyper)
% Calulates statistics of each connected component.
% INPUTS: L - connected components as returned by conncomp_init()
%         L_iso - Same as L but with isotropic voxels for
%                 calculating the compactness and relative anisotropy
%         S_roi - ROI mask (same size as L)
%         S_weight - Weight mask (same size as L)
%         F - Vector that specifies the voxel size
%         SM_hypo - Mask selecting voxels that appear hypointense on T1w volumes
%         SM_hyper - Mask selecting voxels that appear hyperintense on T1w volumes
% OUTPUT: CC - List of structures where each structure contains
%              the label, volume, compactness, location, relative
%              anisotropy, maximal in-plane area, number of slices
%              relative contrast
%
CC = struct( 'lab', -1, ...
             'vol', 0, ...
             'com', NaN, ...
             'loc', -1, ...
             'ra', NaN, ...
             'ma', 0, ...
             'nsl', NaN, ...
             'intvar', NaN, ...
             'phypo', NaN);
idx = 1;
[Lab, Loc] = conncomp_mask(L, S_roi, 0.5, S_weight);
N_lab = length(Lab);
for lab_idx = 1:N_lab
    SM = L==Lab(lab_idx);
    SM_iso = L_iso==Lab(lab_idx);
    if get_volume(SM, F) >= 1
        CC_tmp = struct('lab', Lab(lab_idx), ...
                         'vol', get_volume(SM, F), ...
                         'com', get_compactness(SM_iso), ...
                         'loc', Loc(lab_idx), ...
                         'ra', get_relativeanisotropy(SM_iso), ...
                         'ma', get_max_area(SM, F), ...
                         'nsl', get_nslices(SM), ...
                         'intvar', get_intvar_slice(SM, S_weight), ...
                         'phypo', get_phypo(SM, SM_hypo), ...
                         'phyper', get_phyper(SM, SM_hyper))
        if idx == 1
            CC = CC_tmp;
        else
            CC(idx) = CC_tmp;
        end
        idx = idx + 1;
    end
end


%-------------------------------------------------------------------------
function [phypo] = get_phypo(SM, SM_hypo)
SM_tmp = SM & SM_hypo;
phypo = sum(SM_tmp(:))/sum(SM(:));


%-------------------------------------------------------------------------
function [phyper] = get_phyper(SM, SM_hyper)
SM_tmp = SM & SM_hyper;
phyper = sum(SM_tmp(:))/sum(SM(:));


%-------------------------------------------------------------------------
function [intvar] = get_intvar_slice(SM, S)
cnt = 0;
intvar = [];
for idx_slice = 1:size(SM, 3)
    S_tmp = S(:, :, idx_slice);
    SM_tmp = SM(:, :, idx_slice);
    if sum(SM_tmp(:)) > 0
        cnt = cnt + 1;
        intvar(cnt) = var(double(S_tmp(SM_tmp)));
    end
end
if ~isempty(intvar)
    intvar(intvar == 0) = [];
end
if ~isempty(intvar)
    intvar = median(intvar);
else
    intvar = NaN;
end


%-------------------------------------------------------------------------
function [Nslices] = get_nslices(SM)
Nslices = 0;
for slice = 1:size(SM, 3)
    M_tmp = SM(:, :, slice);
    if sum(M_tmp(:))
        Nslices = Nslices + 1;
    end
end


%-------------------------------------------------------------------------
function [RA] = get_relativeanisotropy(SM)
V = get_volume(SM);
if V > 1
    [X_idx, Y_idx, Z_idx] = ind2sub(size(SM), find(SM));
    I_cov = cov([X_idx, Y_idx, Z_idx]);
    [V, D] = eig(I_cov);
    RA = sqrt(1/2) * ...
         sqrt((D(1,1)-D(2,2))^2 + (D(2,2)-D(3,3))^2 + (D(1,1)-D(3,3))^2) / ...
         (sum(diag(D)));
else
    RA = 0;
end


%-------------------------------------------------------------------------
function [Cd] = get_compactness(SM)
V = get_volume(SM);
if V > 1
    [SM_perim, A] = get_perim_voxels(SM, true);
    Ac = (6*V-A)/2;
    Ac_min = V-1;
    Ac_max = 3*(V - V^(2/3));
    Cd = (Ac - Ac_min)/(Ac_max - Ac_min);
    if Cd < 0
        fprintf('');
    end
else
    Cd = 1;
end


%-------------------------------------------------------------------------
function [Vol, Loc] = get_location_vol(SM, S, F)
Loc = unique(S(:))';
Loc = Loc(2:end);
N_loc = length(Loc);
Vol = zeros(1, N_loc);
S_lab = S .* cast(SM, class(S));
for idx = 1:N_loc
    Vol(idx) = get_volume(S_lab == Loc(idx), F);
end


%-------------------------------------------------------------------------
function [Ma] = get_max_area(SM, F)
Vol = zeros(size(SM, 3), 1);
for slice = 1:size(SM, 3)
    Vol(slice) = get_volume(SM(:, :, slice), [F(1:2) 1]);
end
Ma = max(Vol);
