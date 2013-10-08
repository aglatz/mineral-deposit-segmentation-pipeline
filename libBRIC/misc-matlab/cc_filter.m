function [S_hypos, S_hypos_hypo, S_hypos_hyper] = ...
        cc_filter(S_gre, S_t1w, S_roi, SM_hypos, SM_ntis, intvar_q, I_thr)

L = conncomp_init(SM_hypos, 3);
[Lab, Loc] = conncomp_mask(L, S_roi, 0.5, S_gre);

if ~isempty(intvar_q)
    intvar_thr = get_intvar_ntis(SM_ntis, S_gre, intvar_q);
else
    intvar_thr = NaN;
end

fprintf('Before: %d; intvar_thr=%0.2f\n', sum(SM_hypos(:)), intvar_thr);

N_lab = length(Lab);
S_hypos = zeros(size(S_roi), class(S_roi));
S_hypos_hypo = zeros(size(S_roi), class(S_roi));
S_hypos_hyper = zeros(size(S_roi), class(S_roi));
for lab_idx = 1:N_lab
	SM_cc = L == Lab(lab_idx);
    M = I_thr(:, end) == Loc(lab_idx);
    phypo = get_phypo(SM_cc, S_t1w, I_thr(M, 4));
    phyper = get_phyper(SM_cc, S_t1w, I_thr(M, 5));
    if ~isempty(intvar_q)
        intvar = get_intvar_slice(SM_cc, S_gre);
        res = intvar > intvar_thr;
    else
        intvar = NaN;
        res = true;
    end
    fprintf('Lab:%d loc:%d phypo:%0.2f phyper:%0.2f intvar:%0.2f', ...
        lab_idx, Loc(lab_idx), phypo, phyper, intvar);
    if phypo < 0.5 && phyper < 0.5 && ~(phypo > 0.05 && phyper > 0.05) && res
        % Apply T1w threshold
        SM_cc(SM_cc) = S_t1w(SM_cc) > I_thr(M, 4);
        % all
        S_hypos = S_hypos + cast(SM_cc, class(S_hypos)) .* Loc(lab_idx);
        % hypo hypo
        [~, SM_tmp] = get_phypo(SM_cc, S_t1w, I_thr(M, 2));
        S_hypos_hypo = S_hypos_hypo + cast(SM_tmp, class(S_hypos_hypo)) .* Loc(lab_idx);
        % hypo hyper
        [~, SM_tmp] = get_phyper(SM_cc, S_t1w, I_thr(M, 3));
        S_hypos_hyper = S_hypos_hyper + cast(SM_tmp, class(S_hypos_hyper)) .* Loc(lab_idx);
        fprintf('*\n');
    else
        fprintf('\n');
    end
end

fprintf('After: %d\n', sum(logical(S_hypos(:))));


%-------------------------------------------------------------------------
function [phypo, SM_hypo] = get_phypo(SM, S, I_thr)
SM_hypo = false(size(SM));
SM_hypo(SM) = S(SM) < I_thr;
phypo = sum(SM_hypo(:))/sum(SM(:));

%-------------------------------------------------------------------------
function [phyper, SM_hyper] = get_phyper(SM, S, I_thr)
SM_hyper = false(size(SM));
SM_hyper(SM) = S(SM) > I_thr;
phyper = sum(SM_hyper(:))/sum(SM(:));

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

