function [Ret] = validate_raw(S, S_ref, S_base, F)
% Calculates validation statistics between generated and reference masks:
% - Confusion matrix, and from it
% - Accuracy, Sensitivity, Specificity, Jaccard, Dice, Kappa, Precision,
%   reject Precision
% Additionally, the total volume of the generated and reference masks
% are calculated.
% INPUTS: S - generated mask
%         S_ref - reference mask
%         S_base - base mask. This is needed to calculate
%                  the false positive rate, which is required for e.g.
%                  the sensitivity or kappa calculation. This mask should
%                  be bigger than the generated or reference masks. If
%                  empty this mask will be derived from the generated and
%                  reference masks.
%         F - voxel size
% OUTPUTS: Ret - Structure containing the validation results.
%

% Tests:
% S1=ones(3);S2=ones(3);validate_raw(S1, S2, [], [1 1 1]); => J=1
% S1=ones(3);S2=zeros(3);validate_raw(S1, S2, [], [1 1 1]); => J=0
% S1=eye(3);S2=zeros(3);validate_raw(S1, S2, [], [1 1 1]); => J=0
% S1=eye(3);S2=ones(3);validate_raw(S1, S2, [], [1 1 1]); => J=1/3
% S1=zeros(3);S2=zeros(3);validate_raw(S1, S2, [], [1 1 1]); => J=1

SM = logical(S);
SM_ref = logical(S_ref);

% Create base masks if it does not exist
if ~isempty(S_base)
    SM_base = logical(S_base);
else
	SM_dil = cat(3, [[0 0 0]; [0 1 0]; [0 0 0]], ...
                    [[1 1 1]; [1 1 1]; [1 1 1]], ...
                    [[0 0 0]; [0 1 0]; [0 0 0]]);
	SM_base = imdilate(SM | SM_ref, SM_dil);
end
    
% Everything will be derived from the confusion matrix
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
if (sum(ConfMat(3, 1:2)) - sum(ConfMat(1:2, 3))) > prod(F)
	error('Confusion matrix inconsitency');
end
ConfMat(3, 3) = sum(ConfMat(3, 1:2));


%-------------------------------------------------------------------------
function Ret = get_accuracy(ConfMat)
if ~ConfMat(3, 3)
	Ret = 1; % special case
else
    Ret = (ConfMat(1, 1) + ConfMat(2, 2))/(ConfMat(3, 3));
end


%-------------------------------------------------------------------------
function Ret = get_sensitivity(ConfMat)
Ret = ConfMat(1, 1) / ConfMat(1, 3);


%-------------------------------------------------------------------------
function Ret = get_specificity(ConfMat)
Ret = ConfMat(2, 2) / ConfMat(2, 3);


%-------------------------------------------------------------------------
function Ret = get_jaccard(ConfMat)
if ~ConfMat(3, 3)
    Ret = 1; % special case
else
    Ret = ConfMat(1, 1) / (ConfMat(3, 3) - ConfMat(2, 2));
end


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
if ~ConfMat(3, 3)
    Ret = 1; % special case
else
    po = (ConfMat(1, 1) + ConfMat(2, 2))/ConfMat(3, 3);
    pe = (ConfMat(3, 1)*ConfMat(1, 3))/ConfMat(3, 3).^2 + ...
         (ConfMat(3, 2)*ConfMat(2, 3))/ConfMat(3, 3).^2;
    Ret = (po-pe)/(1-pe);
end

