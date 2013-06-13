function [J, D, V, V_ref, Subjects, FPC, FPV, FPSubjects] = load_matdata(MatName)
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
else
    Subjects_all = Subjects;
end
M_vol = [Ret.Vol] > 0;
M_vol_ref = [Ret.Vol_ref] > 0;

J = [Ret.Jaccard];
J = J(M_vol_ref);

D = [Ret.Dice];
D = D(M_vol_ref);

V_all = [Ret.Vol];
V = V_all(M_vol_ref);

V_ref_all = [Ret.Vol_ref];
V_ref = V_ref_all(M_vol_ref);

Subjects = Subjects_all(M_vol_ref);

FPC = sum(M_vol & ~M_vol_ref);

FPV = V_all(M_vol & ~M_vol_ref);

FPSubjects = Subjects_all(M_vol & ~M_vol_ref);