function [Subjects, J, J0, D, V, V_ref, CNR, q] = load_matdata(MatName)
load(MatName);
if exist('out', 'var')
    Ret = [out.Ret];
    Subjects = cell(length(Ret), 1);
    cnt = 1;
    for idx = 1:length(out)
        Tmp = out(idx).Subjects;
        Subjects(cnt:(cnt+length(Tmp)-1)) = Tmp;
        cnt = cnt + length(Tmp);
    end
end

E = [Ret.edit];
J = [E.Jaccard];
if isfield(Ret, 'thr')
    Tmp = [Ret.thr];
    J0 = [Tmp.Jaccard];
end
D = [E.Dice];
V = [E.Vol];
V_ref = [E.Vol_ref];
if isfield(Ret, 'CNR_oli')
    CNR = [Ret.CNR_oli];
    Input = [Ret.Input];
    q = [Input.intvar_thr];
else
    CNR = NaN(size(J));
    q = NaN(size(J));
end
