function [CC] = conncomp_list(L, L_iso, S_roi, S_weight, F, SM_hypo)
% Calulates statistics of each connected component.
% INPUTS: L - connected components as returned by conncomp_init()
%         L_iso - Same as L but with isotropic voxels for
%                 calculating the compactness and relative anisotropy
%         S_roi - ROI mask (same size as L)
%         S_weight - Weight mask (same size as L)
%         F - Vector that specifies the voxel size
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
             'cont', NaN, ...
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
                         'cont', get_contrast(SM, S_weight), ...
                         'phypo', get_phypo(SM, SM_hypo))
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
function [cont] = get_contrast(SM, S)
hu = 0;
cont = [];
for idx = 1:size(SM, 3)
    S_tmp = S(:, :, idx);
    SM_tmp = SM(:, :, idx);
    if sum(SM_tmp(:)) > 0
        hu = hu + 1;
        Q = quantile(double(S_tmp(SM_tmp)), [.25 .5 .75]);
        cont(hu) = (Q(3)-Q(1)) / Q(2);
    end
end
if ~isempty(cont)
    cont(cont == 0) = [];
end
if ~isempty(cont)
    cont = median(cont);
else
    cont = NaN;
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
function [I_ret] = get_threshold(S, SM, SM_ref, TE)
[I_tmp, I_ref] = volume_stats(S, SM_ref, [], [], []);
deltaR2s = get_deltaR2s(S(SM), I_ref, TE);
I_ret = quantile(deltaR2s, .05); % Selects 95% of the voxels


%-------------------------------------------------------------------------
% function [Ratio] = get_contrast(S_gre, S_t1w, SM, type)
% [I_pm_t1w, I_cm_t1w] = get_3D_pmcm(S_t1w, SM);
% [I_pm_gre, I_cm_gre] = get_3D_pmcm(S_gre, SM);
% switch type
%     case {'cm'},
%         Ratio = I_cm_gre;
%     case {'pm'},
%         Ratio = I_pm_gre;
%     case {'raw'},
%         Ratio = [I_cm_t1w I_pm_t1w I_cm_gre I_pm_gre];
%     otherwise,
%         Ratio = NaN;
% end

%-------------------------------------------------------------------------
function [Res] = get_correlation(I, type)
if size(I, 1) > 2
    try
        switch type
            case {'ang'},
                [tmp, tmp, tmp, tmp, tmp, C] = pcomp_find(I); clear tmp;
                r = C(1,2)/sqrt(C(1, 1)*C(2, 2)); % corrcoef()
                S = std(double(I));
                Res = atand(2*r*S(1)*S(2)/(S(1)^2-S(2)^2))/2;
            case {'cor'},
                Tmp = corr(double(I), 'type', 'Spearman');
                Res = Tmp(1, 2);
            case {'corprob'},
                [Tmp, Tmp] = corr(double(I), 'type', 'Spearman');
                Res = Tmp(1, 2);
            otherwise,
                Res = NaN;
        end
    catch
        Res = NaN;
    end
else
    Res = NaN;
end


%-------------------------------------------------------------------------
function [Mra2D] = get_max_relativeanisotropy_2D(SM)
Ra2D = zeros(size(SM, 3), 1);
for slice = 1:size(SM, 3)
    SM_tmp = SM(:, :, slice);
    if get_volume(SM_tmp) ~= 0
        Ra2D(slice) = get_relativeanisotropy_2D(SM_tmp);
    else
        Ra2D(slice) = -1;
    end
end
Mra2D = max(Ra2D);


%-------------------------------------------------------------------------
function [RA] = get_relativeanisotropy_2D(SM)
V = get_volume(SM);
if V > 1
    [X_idx, Y_idx, Z_idx] = ind2sub(size(SM), find(SM));
    if sum(Z_idx - Z_idx(1)*ones(size(Z_idx))) == 0
        I_cov = cov([X_idx, Y_idx]);
        [V, D] = eig(I_cov);
        RA = abs(D(1,1)-D(2,2))/(D(1,1)+D(2,2));
    else
        RA = -1;
    end
else
    RA = 0;
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
function [Loc] = get_location(SM, S, F)
[Vol, Loc] = get_location_vol(SM, S, F);
if sum(Vol) > 0
    idx = find(Vol == max(Vol));
    Loc = Loc(idx(1));
else
    Loc = -1;
end

%-------------------------------------------------------------------------
function [Ma] = get_max_area(SM, F)
Vol = zeros(size(SM, 3), 1);
for slice = 1:size(SM, 3)
    Vol(slice) = get_volume(SM(:, :, slice), [F(1:2) 1]);
end
Ma = max(Vol);
