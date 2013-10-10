function [Ret] = get_tis_ints(S, SM_ntis, SM_feat, SM_valid)
% Normal appearing tissue
SM_tmp1 = SM_ntis & SM_valid;
Ret.I_ntis = quantile(S(SM_tmp1), [.05 .25 .5 .75 .95]);
Ret.I_ntis_fvalid = sum(SM_tmp1(:))/sum(SM_ntis(:));
% Features
SM_tmp2 = SM_feat & SM_valid;
Ret.I_feat = quantile(S(SM_tmp2), [.05 .25 .5 .75 .95]);
Ret.I_feat_fvalid = sum(SM_tmp2(:))/sum(SM_feat(:));
% Appearance
a = 0.05;
Ret.P_feat = ranksum(S(SM_tmp1), S(SM_tmp2), 'alpha', a);
Ret.A_feat = get_appearance(Ret.I_ntis(3), Ret.I_feat(3), Ret.P_feat, a);