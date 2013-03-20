function [Ret] = validate(path, name, name_ref, name_base, name_voi, use_voi)
% Calculates validation statistics between generated and reference masks:
% - Confusion matrix, and from it
% - Accuracy, Sensitivity, Specificity, Jaccard, Dice, Kappa, Precision,
%   reject Precision
% Additionally, the total volume of the generated and reference masks
% are calculated.
% INPUTS: path - subject path
%         name - filename of the generated mask (relative to subject dir)
%         name_ref - filename of the reference mask (relative to subject dir)
%         name_base - filename of the base mask. This is needed to calculate
%                     the false positive rate, which is required for e.g.
%                     the sensitivity or kappa calculation. This mask should
%                     be bigger than the generated or reference masks. If
%                     empty this mask will be derived from the generated and
%                     reference masks.
%         name_voi - filename of the VOI masks to speed up processing by
%                    ignoring blank slices (see roi_init()).
%         use_voi - flag indicating whether name_voi should be used or not.
% OUTPUTS: Ret - Structure containing the validation results.
%

if use_voi
    S_voi = load_series([path '/' name_voi], []);
    Roi = roi_init(S_voi);
	slices = roi_nifti_sliceno(Roi, []);   
else
	slices = [];
end

% New mask
SM = logical(load_series([path '/' name], slices));
NII = load_series([path '/' name], 0);
F = NII.hdr.dime.pixdim(2:4);

% Total gold standard mask
SM_ref = logical(load_series([path '/' name_ref], slices));

if ~isempty(name_base)
	% Approx. of gold standard mask (always bigger than gold
    % standard, but smaller than ROI).
	SM_base = logical(load_series([path '/' name_base], slices));
else
    SM_dil = cat(3, [[0 0 0]; [0 1 0]; [0 0 0]], ...
                    [[1 1 1]; [1 1 1]; [1 1 1]], ...
                    [[0 0 0]; [0 1 0]; [0 0 0]]);
	SM_base = imdilate(SM | SM_ref, SM_dil);
end

% Confusion matrix
ConfMat = get_confmat(SM, SM_ref, SM_base, F);

Ret = struct;
Ret.Accuracy = get_accuracy(ConfMat);
Ret.Sensitivity = get_sensitivity(ConfMat);
Ret.Specificity = get_specificity(ConfMat);
Ret.Jaccard = get_jaccard(ConfMat);
Ret.Dice = get_dice(ConfMat);
Ret.Kappa = get_kappa(ConfMat);
Ret.Precision = get_precision(ConfMat);
Ret.RejPrecision = get_rejprecision(ConfMat);
Ret.Vol = get_volume(SM, F);
Ret.Vol_ref = get_volume(SM_ref, F);
Ret.Vol_base = get_volume(SM&SM_base, F);
Ret.Vol_ref_base = get_volume(SM_ref&SM_base, F);


fprintf(1, '\n------ Classification statistics ------\n');
ConfMat
Ret
fprintf(1, '-----------------------------------------\n');


%-------------------------------------------------------------------------
function [ConfMat] = get_confmat(SM, SM_ref, SM_base, F)
ConfMat = zeros(3, 3);
ConfMat(1, 1) = get_volume(SM & SM_ref, F); %TP
ConfMat(2, 2) = get_volume((~SM & SM_base) & (~SM_ref & SM_base), F); %TN
ConfMat(1, 2) = get_volume((SM_ref & SM_base) & ~(SM & SM_ref), F); %FN
ConfMat(2, 1) = get_volume((SM & SM_base) & ~(SM & SM_ref), F); %FP
ConfMat(3, 1:2) = sum(ConfMat(1:2, 1:2), 1);
ConfMat(1:2, 3) = sum(ConfMat(1:2, 1:2), 2);
if sum(ConfMat(3, 1:2)) ~= sum(ConfMat(1:2, 3))
	error('Confusion matrix inconsitency');
end
ConfMat(3, 3) = sum(ConfMat(3, 1:2));


%-------------------------------------------------------------------------
function Ret = get_accuracy(ConfMat)
Ret = (ConfMat(1, 1) + ConfMat(2, 2))/(ConfMat(3, 3));


%-------------------------------------------------------------------------
function Ret = get_sensitivity(ConfMat)
Ret = ConfMat(1, 1) / ConfMat(1, 3);


%-------------------------------------------------------------------------
function Ret = get_specificity(ConfMat)
Ret = ConfMat(2, 2) / ConfMat(2, 3);


%-------------------------------------------------------------------------
function Ret = get_jaccard(ConfMat)
Ret = ConfMat(1, 1) / (ConfMat(3, 3) - ConfMat(2, 2));


%-------------------------------------------------------------------------
function Ret = get_precision(ConfMat)
Ret = ConfMat(1, 1) / ConfMat(3, 1);


%-------------------------------------------------------------------------
function Ret = get_rejprecision(ConfMat)
Ret = ConfMat(2, 2) / ConfMat(3, 2);

%-------------------------------------------------------------------------
function Ret = get_dice(ConfMat)
Ret = 2*get_jaccard(ConfMat)/(get_jaccard(ConfMat)+1);

%-------------------------------------------------------------------------
function Ret = get_kappa(ConfMat)
po = (ConfMat(1, 1) + ConfMat(2, 2))/ConfMat(3, 3);
pe = (ConfMat(3, 1)*ConfMat(1, 3))/ConfMat(3, 3).^2 + ...
     (ConfMat(3, 2)*ConfMat(2, 3))/ConfMat(3, 3).^2;
Ret = (po-pe)/(1-pe);

