function [J, D, V, V_ref, Subjects, FPC, FPV, FPSubjects, Sens, Spec, Conf] = load_matdata(MatName, tpflag)
load(MatName);
if exist('out', 'var')
    Ret = [out.Ret];
    Subjects_all = cell(length(Ret), 1);
    cnt = 1;
    for idx = 1:length(out)
        Tmp = out(idx).Subjects;
        Subjects_all(cnt:(cnt+length(Tmp)-1)) = Tmp;
        cnt = cnt + length(Tmp);
    end
%    P = reshape([out.P], 2, length(out))';
else
    Subjects_all = Subjects;
%    P = [1 0];
end
M_vol = [Ret.Vol] > 0;
M_vol_ref = [Ret.Vol_ref] > 0;

J = [Ret.Jaccard];
if tpflag
	J = J(M_vol & M_vol_ref);
end

D = [Ret.Dice];
if tpflag
	D = D(M_vol & M_vol_ref);
end

V_all = [Ret.Vol];
if tpflag
    V = V_all(M_vol & M_vol_ref);
else
    V = V_all;
end

V_ref_all = [Ret.Vol_ref];
if tpflag
    V_ref = V_ref_all(M_vol & M_vol_ref);
else
    V_ref = V_ref_all;
end

Subjects = Subjects_all(M_vol & M_vol_ref);

FPC = sum(M_vol & ~M_vol_ref);

FPV = V_all(M_vol & ~M_vol_ref);

FPSubjects = Subjects_all(M_vol & ~M_vol_ref);

% Conf(1,1) = TN, Conf(2,2) = TP, Conf(1,2) = FP, Conf(2,1) = FN
Conf = confusionmat(M_vol_ref, M_vol);
Sens = Conf(2,2)/sum(Conf(2,:));
Spec = Conf(1,1)/sum(Conf(1,:));